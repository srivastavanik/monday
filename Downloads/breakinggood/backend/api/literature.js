const express = require("express");
const axios = require("axios");
const xml2js = require("xml2js");
const router = express.Router();
const logger = require("../utils/logger"); // Added logger

// PubMed API base URLs
const PUBMED_ESEARCH_URL =
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
const PUBMED_EFETCH_URL =
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const PUBMED_ELINK_URL =
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi";
const BIOC_PMC_URL =
  "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json";

// Crossref API base URL
const CROSSREF_API_URL = "https://api.crossref.org/works";

// Search literature by keyword
router.get("/search", async (req, res) => {
  try {
    const { query, limit = 10, page = 1 } = req.query;

    if (!query) {
      return res.status(400).json({ error: "Search query is required" });
    }

    logger.info(
      `Searching Crossref for: ${query}, limit: ${limit}, page: ${page}`
    );

    // Build the query with valid filters and selects
    const searchParams = {
      query: query,
      rows: limit,
      offset: (page - 1) * limit,
      mailto: "your-email@example.com", // Required by Crossref API
      filter: [
        "has-abstract:true", // Only include articles with abstracts
        "type:journal-article", // Focus on journal articles
        "has-references:true", // Only include articles with references
        "has-license:true", // Only include articles with licenses
        "has-affiliation:true", // Only include articles with author affiliations
      ].join(","),
      select: [
        "DOI",
        "title",
        "abstract",
        "author",
        "container-title",
        "created",
        "URL",
        "type",
        "reference",
        "is-referenced-by-count",
        "license",
        "funder",
        "publisher",
        "subject",
        "published",
        "ISSN",
      ].join(","),
    };

    const response = await axios.get(CROSSREF_API_URL, {
      params: searchParams,
    });

    const items = response.data.message.items || [];
    const totalResults = response.data.message["total-results"];

    const formattedArticles = items
      .map((article) => {
        // Verify author information
        const authors =
          article.author
            ?.map((author) => {
              const name = `${author.family || ""}, ${
                author.given || ""
              }`.trim();
              const affiliation = author.affiliation?.[0]?.name || "";
              return { name, affiliation };
            })
            .filter((author) => author.name) || [];

        // Get funding information
        const funders =
          article.funder?.map((f) => ({
            name: f.name,
            award: f.award?.[0],
          })) || [];

        // Get publisher information
        const publisher = article.publisher || "N/A";

        // Get subjects/categories
        const subjects = article.subject || [];

        return {
          doi: article.DOI,
          title: article.title?.[0] || "No Title Available",
          abstract: article.abstract || "",
          authors: authors
            .map(
              (a) => `${a.name}${a.affiliation ? ` (${a.affiliation})` : ""}`
            )
            .join(", "),
          journal: article["container-title"]?.[0] || "N/A",
          publicationYear: article.created?.["date-parts"]?.[0]?.[0] || "N/A",
          url: article.URL,
          type: article.type,
          references: article.reference?.length || 0,
          citations: article["is-referenced-by-count"] || 0,
          license: article.license?.[0]?.URL || null,
          funders: funders,
          publisher: publisher,
          subjects: subjects,
          metadata: {
            hasAbstract: !!article.abstract,
            hasAuthors: authors.length > 0,
            hasReferences: article.reference?.length > 0,
            hasLicense: !!article.license?.[0],
            hasFunding: funders.length > 0,
            hasPublisher: !!publisher,
            hasSubjects: subjects.length > 0,
          },
        };
      })
      .filter(
        (article) =>
          // Additional quality checks
          article.metadata.hasAbstract &&
          article.metadata.hasAuthors &&
          article.metadata.hasReferences &&
          article.authors.length > 0
      );

    logger.info(
      `Found ${formattedArticles.length} articles (total ${totalResults})`
    );

    return res.json({
      results: formattedArticles,
      total: totalResults,
      page: parseInt(page),
      limit: parseInt(limit),
      filters: {
        hasAbstract: true,
        hasReferences: true,
        hasLicense: true,
        hasAffiliation: true,
        type: "journal-article",
      },
    });
  } catch (error) {
    logger.error(`Literature search error: ${error.message}`, {
      stack: error.stack,
      response: error.response?.data,
    });
    return res.status(500).json({
      error: "Error searching literature",
      details: error.message,
    });
  }
});

// Get article details by PMID or PMCID, including full text if available via BioC
router.get("/pubmed/:id", async (req, res) => {
  // Changed endpoint to /pubmed/:id
  try {
    const { id } = req.params;
    const isPmcId = id.toUpperCase().startsWith("PMC");

    if (!id) {
      return res.status(400).json({ error: "PubMed ID or PMC ID is required" });
    }

    logger.info(`Fetching details for ID: ${id}`);

    // 1. Fetch standard details from eFetch first
    let articleDetails = {};
    try {
      const fetchResponse = await axios.get(PUBMED_EFETCH_URL, {
        params: {
          db: "pubmed",
          id: id,
          retmode: "xml",
        },
      });

      const parser = new xml2js.Parser({
        explicitArray: false,
        ignoreAttrs: true,
      });
      const result = await parser.parseStringPromise(fetchResponse.data);

      if (!result.PubmedArticleSet?.PubmedArticle) {
        throw new Error("Article not found via eFetch");
      }

      const article = result.PubmedArticleSet.PubmedArticle;
      const medlineCitation = article.MedlineCitation;
      const articleData = medlineCitation.Article;
      const pubmedData = article.PubmedData;

      let authors = [];
      if (articleData.AuthorList && articleData.AuthorList.Author) {
        const authorList = Array.isArray(articleData.AuthorList.Author)
          ? articleData.AuthorList.Author
          : [articleData.AuthorList.Author];
        authors = authorList
          .map((author) =>
            author.LastName && author.ForeName
              ? `${author.LastName}, ${author.ForeName}`
              : author.LastName || author.CollectiveName || ""
          )
          .filter(Boolean);
      }

      let abstract = "";
      if (articleData.Abstract && articleData.Abstract.AbstractText) {
        abstract = Array.isArray(articleData.Abstract.AbstractText)
          ? articleData.Abstract.AbstractText.map((t) =>
              typeof t === "string" ? t : t._
            ).join(" ")
          : typeof articleData.Abstract.AbstractText === "string"
          ? articleData.Abstract.AbstractText
          : articleData.Abstract.AbstractText._ || "";
      }

      const pmcId = pubmedData?.ArticleIdList?.ArticleId?.find(
        (id) => id.IdType === "pmc"
      )?._;
      const doi = pubmedData?.ArticleIdList?.ArticleId?.find(
        (id) => id.IdType === "doi"
      )?._;
      const pmidActual = medlineCitation.PMID;

      articleDetails = {
        pmid: pmidActual,
        pmcid: pmcId,
        title: articleData.ArticleTitle || "No Title Available",
        abstract,
        authors: authors.join(", "),
        journal: articleData.Journal?.Title || "N/A",
        publicationYear:
          articleData.Journal?.JournalIssue?.PubDate?.Year || "N/A",
        doi: doi,
        keywords: medlineCitation.KeywordList?.Keyword
          ? (Array.isArray(medlineCitation.KeywordList.Keyword)
              ? medlineCitation.KeywordList.Keyword
              : [medlineCitation.KeywordList.Keyword]
            ).map((k) => (typeof k === "string" ? k : k._))
          : [],
        meshTerms: medlineCitation.MeshHeadingList?.MeshHeading
          ? (Array.isArray(medlineCitation.MeshHeadingList.MeshHeading)
              ? medlineCitation.MeshHeadingList.MeshHeading
              : [medlineCitation.MeshHeadingList.MeshHeading]
            ).map((term) => term.DescriptorName?._ || term.DescriptorName)
          : [],
      };
    } catch (efetchError) {
      logger.error(`eFetch failed for ${id}: ${efetchError.message}`);
      // Continue to try BioC if possible
      if (!isPmcId && !articleDetails.pmcid) {
        // If we don't have a PMCID, we can't use BioC easily, return error
        return res.status(404).json({
          error: "Article details not found via eFetch and PMCID is unknown.",
        });
      }
    }

    // 2. Attempt to fetch full text from BioC using PMCID if available
    const idForBioC = isPmcId ? id : articleDetails.pmcid;
    let fullText = null;
    let sections = [];
    if (idForBioC) {
      try {
        logger.info(`Attempting to fetch BioC JSON for ID: ${idForBioC}`);
        const biocResponse = await axios.get(
          `${BIOC_PMC_URL}/${idForBioC}/unicode`
        );
        const biocData = biocResponse.data;

        // Extract passages (text sections)
        if (biocData && biocData.documents && biocData.documents.length > 0) {
          const document = biocData.documents[0];
          fullText = ""; // Concatenate all text for simple fullText field
          document.passages?.forEach((passage) => {
            const sectionTitle =
              passage.infons?.section_type || "Unknown Section";
            const text = passage.text || "";
            sections.push({ title: sectionTitle, text: text });
            fullText += text + "\n\n";
          });
          logger.info(
            `Successfully fetched and parsed BioC data for ${idForBioC}`
          );
        } else {
          logger.warn(
            `BioC data structure unexpected or empty for ${idForBioC}`
          );
        }
      } catch (biocError) {
        if (biocError.response?.status === 404) {
          logger.warn(`BioC full text not found for ID: ${idForBioC}`);
        } else {
          logger.error(
            `Error fetching BioC full text for ${idForBioC}: ${biocError.message}`
          );
        }
        // Full text is optional, so we continue without it
      }
    }

    // 3. Fetch related articles (remains the same)
    let relatedPmids = [];
    try {
      const relatedResponse = await axios.get(PUBMED_ELINK_URL, {
        params: {
          dbfrom: "pubmed",
          db: "pubmed",
          id: articleDetails.pmid || id, // Use fetched PMID if available
          cmd: "neighbor_score",
          retmode: "json",
        },
      });
      const linkSetDb = relatedResponse.data.linksets?.[0]?.linksetdbs?.find(
        (db) => db.linkname === "pubmed_pubmed"
      );
      relatedPmids = linkSetDb ? linkSetDb.links?.slice(0, 5) || [] : [];
    } catch (e) {
      logger.error(`Error parsing related articles for ${id}: ${e.message}`);
    }

    // Combine all fetched data
    const finalDetails = {
      ...articleDetails,
      fullText: fullText?.trim() || null,
      sections: sections.length > 0 ? sections : null,
      relatedArticles: relatedPmids,
    };

    // If we couldn't get basic details from eFetch but got BioC, fill from BioC if possible
    if (!articleDetails.title && fullText) {
      // Attempt basic extraction from BioC structure if needed (less reliable)
      const biocDoc = biocResponse?.data?.documents?.[0];
      articleDetails.title =
        biocDoc?.passages?.find((p) => p.infons?.type === "title")?.text ||
        "Title Not Found";
      // Add other extractions if necessary
    }

    return res.json(finalDetails);
  } catch (error) {
    logger.error(
      `Error fetching article details for ${req.params.id}: ${error.message}`,
      { stack: error.stack, response: error.response?.data }
    );
    return res.status(500).json({
      error: "Error fetching article details",
      details: error.message,
    });
  }
});

// Analyze literature using Claude
router.post("/analyze", async (req, res) => {
  try {
    const { articles, query } = req.body; // articles can be PMIDs, abstracts, or full texts

    if (
      !articles ||
      !Array.isArray(articles) ||
      articles.length === 0 ||
      !query
    ) {
      return res
        .status(400)
        .json({ error: "Articles (list) and analysis query are required." });
    }

    // Prepare context for Claude
    let context = "Analyze the following literature excerpts:\n\n";
    articles.forEach((article, index) => {
      context += `--- Article ${index + 1} ---\n`;
      if (typeof article === "string") {
        // Assume it's abstract or text
        context += `${article.substring(0, 1500)}...\n`; // Limit context size
      } else if (article.abstract) {
        context += `Title: ${
          article.title
        }\nAbstract: ${article.abstract.substring(0, 1500)}...\n`;
      } else if (article.pmid) {
        context += `PMID: ${article.pmid}\nTitle: ${article.title || "N/A"}\n`;
      }
      context += `\n`;
    });

    const userPrompt = `${query}\n\nBased on the provided literature context:\n${context}`;

    // Call Claude via AI service
    const analysisResponse = await axios.post(
      "http://localhost:5000/api/ai/ask",
      {
        question: query, // Keep the original query separate for clarity
        context: context, // Provide the compiled context
      }
    );

    res.json({ analysis: analysisResponse.data.response });
  } catch (error) {
    logger.error(`Error analyzing literature: ${error.message}`, {
      stack: error.stack,
    });
    res.status(500).json({
      error: "Error analyzing literature",
      details: error.message,
    });
  }
});

// Create research note
router.post("/notes", async (req, res) => {
  try {
    const { articleId, content, title } = req.body;

    if (!articleId || !content) {
      return res
        .status(400)
        .json({ error: "Article ID and content are required" });
    }

    // Here you would typically save the note to your database
    // For now, we'll just return a success response
    const note = {
      id: Date.now().toString(),
      articleId,
      title: title || "Untitled Note",
      content,
      createdAt: new Date().toISOString(),
    };

    return res.status(201).json(note);
  } catch (error) {
    logger.error(`Error creating research note: ${error.message}`);
    return res.status(500).json({
      error: "Error creating research note",
      details: error.message,
    });
  }
});

module.exports = router;

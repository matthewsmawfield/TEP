#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

/**
 * HTML to Markdown Converter for the TEP theory site.
 * Reads manifest.json + components/*.html directly (no dist dependency)
 * and writes 0-TEP-v{version}-{codename}.md at the repo root.
 */
class HTMLToMarkdownConverter {
    decodeEntities(text) {
        return text
            .replace(/&amp;/g, '&')
            .replace(/&lt;/g, '<')
            .replace(/&gt;/g, '>')
            .replace(/&quot;/g, '"')
            .replace(/&#39;/g, "'")
            .replace(/&nbsp;/g, ' ')
            .replace(/&mdash;/g, '—')
            .replace(/&ndash;/g, '–')
            .replace(/&sect;/g, '§')
            .replace(/&ldquo;/g, '“')
            .replace(/&rdquo;/g, '”')
            .replace(/&lsquo;/g, '‘')
            .replace(/&rsquo;/g, '’')
            .replace(/&hellip;/g, '…')
            .replace(/&times;/g, '×')
            .replace(/&minus;/g, '−')
            .replace(/&plusmn;/g, '±')
            .replace(/&deg;/g, '°')
            .replace(/&thinsp;/g, ' ')
            .replace(/&sup2;/g, '²')
            .replace(/&sup3;/g, '³')
            .replace(/&sup1;/g, '¹')
            .replace(/&#10004;/g, '✔')
            .replace(/&#10008;/g, '✘');
    }

    cleanInlineHtml(html) {
        return this.decodeEntities(
            html
                .replace(/<(strong|b)\b[^>]*>([\s\S]*?)<\/(strong|b)\s*>/gi, '**$2**')
                .replace(/<(em|i)\b[^>]*>([\s\S]*?)<\/(em|i)\s*>/gi, '*$2*')
                .replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`')
                .replace(/<br\s*\/?>/gi, ' ')
                .replace(/<sub[^>]*>([\s\S]*?)<\/sub>/gi, '<sub>$1</sub>')
                .replace(/<sup[^>]*>([\s\S]*?)<\/sup>/gi, '<sup>$1</sup>')
                .replace(/<\/?(?!sub|sup)[a-zA-Z][^>]*>/g, '')
        )
            .replace(/\s+/g, ' ')
            .trim();
    }

    tableToMarkdown(tableHtml) {
        const captionMatch = tableHtml.match(/<caption[^>]*>([\s\S]*?)<\/caption>/i);
        const caption = captionMatch ? this.cleanInlineHtml(captionMatch[1]) : '';
        const rows = [];
        const rowRegex = /<tr[^>]*>([\s\S]*?)<\/tr>/gi;
        let rowMatch;

        while ((rowMatch = rowRegex.exec(tableHtml)) !== null) {
            const cells = [];
            const cellRegex = /<t[hd][^>]*>([\s\S]*?)<\/t[hd]>/gi;
            let cellMatch;

            while ((cellMatch = cellRegex.exec(rowMatch[1])) !== null) {
                cells.push(this.cleanInlineHtml(cellMatch[1]).replace(/\|/g, '\\|'));
            }

            if (cells.length) rows.push(cells);
        }

        if (!rows.length) return '';

        const header = rows[0];
        const normalizeRow = (row) => `| ${header.map((_, index) => row[index] || '').join(' | ')} |`;

        let markdown = '';
        if (caption) markdown += `\n\n${caption}\n\n`;
        markdown += `${normalizeRow(header)}\n`;
        markdown += `| ${header.map(() => '---').join(' | ')} |\n`;
        for (const row of rows.slice(1)) markdown += `${normalizeRow(row)}\n`;
        return `\n\n${markdown}\n`;
    }

    htmlToMarkdown(html) {
        html = html.replace(/<script\b[^<]*(?:(?!<\/script>)<[^<]*)*<\/script>/gis, '');
        html = html.replace(/<style\b[^<]*(?:(?!<\/style>)<[^<]*)*<\/style>/gis, '');
        html = html.replace(/<!--[\s\S]*?-->/g, '');
        html = html.replace(/<nav\b[\s\S]*?<\/nav>/gi, '');

        // Protect math before any tag processing: LaTeX may contain '<'/'|'
        // that would otherwise be misread as markup or table separators.
        // Display blocks are normalised to a single line.
        const mathBlocks = [];
        html = html.replace(/\$\$[\s\S]*?\$\$|\$[^$\n]+?\$/g, (match) => {
            let math = match;
            if (math.startsWith('$$')) {
                const inner = math.slice(2, -2).replace(/\s+/g, ' ').trim();
                math = `$$${inner}$$`;
            }
            mathBlocks.push(math);
            return `@@MATH_${mathBlocks.length - 1}@@`;
        });

        html = html.replace(/<pre[^>]*>\s*<code(?:\s+class=["']language-([^"']+)["'])?[^>]*>([\s\S]*?)<\/code>\s*<\/pre>/gi, (match, lang, code) => {
            const language = (lang || '').trim();
            const decodedCode = this.decodeEntities(code).replace(/\n+$/g, '');
            return `\n\n@@@CODEBLOCK_START:${language}@@@\n${decodedCode}\n@@@CODEBLOCK_END@@@\n\n`;
        });
        html = html.replace(/<pre[^>]*>([\s\S]*?)<\/pre>/gi, (match, code) => {
            const decodedCode = this.decodeEntities(code.replace(/<[^>]+>/g, '')).replace(/\n+$/g, '');
            return `\n\n@@@CODEBLOCK_START:@@@\n${decodedCode}\n@@@CODEBLOCK_END@@@\n\n`;
        });

        html = html.replace(/<table[^>]*>[\s\S]*?<\/table>/gi, (match) => this.tableToMarkdown(match));

        html = html.replace(/<h1[^>]*>([\s\S]*?)<\/h1>/gi, '\n# $1\n\n');
        html = html.replace(/<h2[^>]*>([\s\S]*?)<\/h2>/gi, '\n## $1\n');
        html = html.replace(/<h3[^>]*>([\s\S]*?)<\/h3>/gi, '\n### $1\n');
        html = html.replace(/<h4[^>]*>([\s\S]*?)<\/h4>/gi, '\n#### $1\n\n');
        html = html.replace(/<h5[^>]*>([\s\S]*?)<\/h5>/gi, '\n##### $1\n\n');
        html = html.replace(/<h6[^>]*>([\s\S]*?)<\/h6>/gi, '\n###### $1\n\n');

        html = html.replace(/<figcaption[^>]*>([\s\S]*?)<\/figcaption>/gi, '\n\n$1\n');
        html = html.replace(/<figure[^>]*>([\s\S]*?)<\/figure>/gi, '\n\n$1\n\n');

        html = html.replace(/<img[^>]+>/gi, (tag) => {
            const src = (tag.match(/src=["']([^"']+)["']/i) || [])[1] || '';
            const alt = (tag.match(/alt=["']([^"']*)["']/i) || [])[1] || '';
            // Site-relative paths resolve under site/ when the markdown is
            // rendered from the repository root (e.g. on GitHub).
            const mdSrc = src.replace(/^figures\//, 'site/figures/');
            return src ? `\n![${alt}](${mdSrc})\n` : '';
        });

        html = html.replace(/<blockquote[^>]*>([\s\S]*?)<\/blockquote>/gi, '\n> $1\n\n');
        html = html.replace(/<p[^>]*>([\s\S]*?)<\/p>/gi, (match, content) => {
            const stripped = content.split('\n').map((line) => line.trim()).filter((line) => line.length > 0).join(' ').replace(/ {2,}/g, ' ').trim();
            return `\n\n${stripped}\n\n`;
        });

        html = html.replace(/<(strong|b)\b[^>]*>([\s\S]*?)<\/(strong|b)\s*>/gi, '**$2**');
        html = html.replace(/<(em|i)\b[^>]*>([\s\S]*?)<\/(em|i)\s*>/gi, '*$2*');
        html = html.replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`');

        html = html.replace(/<li[^>]*>\s*([\s\S]*?)\s*<\/li>/gi, '- $1\n');
        html = html.replace(/<br\s*\/?>/gi, '\n');
        html = html.replace(/<hr\s*\/?>/gi, '\n---\n');

        html = html.replace(/<sub[^>]*>([\s\S]*?)<\/sub>/gi, '<sub>$1</sub>');
        html = html.replace(/<sup[^>]*>([\s\S]*?)<\/sup>/gi, '<sup>$1</sup>');
        html = html.replace(/<\/?(?!sub|sup)[a-zA-Z][^>]*>/g, '');
        html = this.decodeEntities(html);

        // Restore code blocks
        html = html.replace(/@@@CODEBLOCK_START:([^@]*)@@@\s*([\s\S]*?)\s*@@@CODEBLOCK_END@@@/g, (match, lang, code) => {
            const language = lang.trim();
            return `\n\n\`\`\`${language}\n${code}\n\`\`\`\n\n`;
        });

        // Restore math (decoding entities that were protected inside $...$)
        html = html.replace(/@@MATH_(\d+)@@/g, (match, idx) => this.decodeEntities(mathBlocks[Number(idx)]));

        // Final cleanup: strip leading spaces, collapse blank lines
        html = html.split('\n').map((line) => line.replace(/^\s+/, '')).join('\n');
        return html.replace(/\n{3,}/g, '\n\n').trim();
    }

    /**
     * Get version info from VERSION.json for filename generation
     */
    getVersionInfo() {
        try {
            const versionPath = path.join(__dirname, '..', 'VERSION.json');
            if (fs.existsSync(versionPath)) {
                const versionData = JSON.parse(fs.readFileSync(versionPath, 'utf8'));
                return {
                    version: versionData.version || '0.0',
                    codename: versionData.codename || 'Unknown',
                };
            }
        } catch (e) {
            console.warn('⚠️  Could not read VERSION.json, using defaults');
        }
        return { version: '0.0', codename: 'Unknown' };
    }

    generateMarkdownFilename() {
        const { version, codename } = this.getVersionInfo();
        return `0-TEP-v${version}-${codename}.md`;
    }

    async convertSiteToMarkdown() {
        console.log('🔄 Converting components to markdown (TEP)...');

        const manifestPath = path.join(__dirname, 'manifest.json');
        if (!fs.existsSync(manifestPath)) throw new Error('manifest.json not found.');
        const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
        const sections = manifest.sections.slice().sort((a, b) => a.order - b.order);

        let allHtml = '';
        for (const section of sections) {
            const componentPath = path.join(__dirname, 'components', section.file);
            if (fs.existsSync(componentPath)) {
                const html = fs.readFileSync(componentPath, 'utf8');
                allHtml += `\n<!-- SECTION: ${section.title} -->\n${html}\n`;
                console.log(`  ✓ ${section.file} (${(html.length / 1024).toFixed(1)} KB)`);
            } else {
                console.warn(`  ⚠ Missing: ${section.file}`);
            }
        }

        console.log(`  Total HTML: ${(allHtml.length / 1024).toFixed(1)} KB`);
        const markdown = this.htmlToMarkdown(allHtml);

        const title = manifest.title || 'Untitled';
        const author = manifest.author || 'Matthew Lukin Smawfield';
        const version = manifest.version || '';
        const firstPublished = manifest.first_published || '';
        const lastUpdated = manifest.date || '';
        const doi = manifest.doi || '';

        const header = `# ${title}
**${author}**
Version: ${version}
First published: ${firstPublished} · Last updated: ${lastUpdated}
DOI: ${doi}

---

`;
        const footer = `
---

*This document was automatically generated from the TEP research site. For the interactive version with figures and enhanced formatting, visit: https://mlsmawfield.com/tep/theory/*

*Source code and data available at: https://github.com/matthewsmawfield/TEP*
`;
        const finalMarkdown = header + markdown + '\n' + footer;

        const filename = this.generateMarkdownFilename();
        const outputPath = path.join(__dirname, '..', filename);
        fs.writeFileSync(outputPath, finalMarkdown, 'utf8');
        console.log(`✅ Markdown saved to: ${outputPath} (${(finalMarkdown.length / 1024).toFixed(1)} KB)`);

        // Keep the shared collection archive (../manuscripts/) in sync
        try {
            const sharedArchiveDir = path.join(__dirname, '..', '..', 'manuscripts');
            if (fs.existsSync(sharedArchiveDir)) {
                const archiveName = path.basename(outputPath);
                const paperPrefix = archiveName.split('-')[0];
                for (const staleFile of fs.readdirSync(sharedArchiveDir)) {
                    if (staleFile !== archiveName && staleFile.endsWith('.md') && staleFile.startsWith(`${paperPrefix}-TEP`)) {
                        fs.rmSync(path.join(sharedArchiveDir, staleFile));
                    }
                }
                fs.copyFileSync(outputPath, path.join(sharedArchiveDir, archiveName));
                console.log(`📄 Copied to shared archive: manuscripts/${archiveName}`);
            }
        } catch (archiveError) {
            console.warn(`⚠️  Could not update shared manuscripts archive: ${archiveError.message}`);
        }
        return outputPath;
    }
}

if (require.main === module) {
    new HTMLToMarkdownConverter().convertSiteToMarkdown().catch((error) => {
        console.error('❌ Markdown conversion failed:', error.message);
        process.exit(1);
    });
}

module.exports = { HTMLToMarkdownConverter };

#!/usr/bin/env python3
"""
Generate PDF from TEP Main Paper
=================================

Generates a high-quality PDF from the main TEP theory paper index.html.
Uses the html_to_pdf.py converter with optimized settings for academic manuscripts.

Usage:
    python scripts/generate_site_pdf.py
    python scripts/generate_site_pdf.py --quality maximum --wait-time 5
"""

import asyncio
import argparse
import subprocess
import sys
import re
from pathlib import Path

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent / 'utils'))

try:
    import yaml
except ImportError:
    yaml = None

from html_to_pdf import HTMLToPDFConverter, create_preset_configs


def load_citation_metadata():
    """Load version and codename from CITATION.cff."""
    base_dir = Path(__file__).parent.parent
    citation_file = base_dir / 'CITATION.cff'
    
    if not citation_file.exists():
        print("⚠️  CITATION.cff not found, using defaults")
        return {'version': '0.9', 'codename': 'Jakarta', 'title': 'Temporal Equivalence Principle'}
    
    try:
        if yaml:
            with open(citation_file, 'r') as f:
                data = yaml.safe_load(f)
            version_str = data.get('version', 'v0.9')
        else:
            # Parse manually if yaml not available
            with open(citation_file, 'r') as f:
                content = f.read()
            version_match = re.search(r'version:\s*"?([^"\n]+)"?', content)
            version_str = version_match.group(1).strip() if version_match else 'v0.9'
        
        # Parse version string like 'v0.9 (Jakarta)'
        pattern = r'^(v?[\d.]+)(?:\s*\(([^)]+)\))?$'
        match = re.match(pattern, version_str.strip())
        
        if match:
            version = match.group(1).lstrip('v')
            codename = match.group(2) or 'Jakarta'
        else:
            version = version_str.lstrip('v')
            codename = 'Jakarta'
        
        return {'version': version, 'codename': codename}
        
    except Exception as e:
        print(f"⚠️  Error parsing CITATION.cff: {e}, using defaults")
        return {'version': '0.9', 'codename': 'Jakarta', 'title': 'Temporal Equivalence Principle'}


def strip_trailing_blank_pages(pdf_path: Path):
    """Remove trailing pages that contain only the page-number footer.

    Chromium printToPDF can emit a final blank page when the body height
    marginally overflows. A page is treated as blank only if its extracted
    text reduces to the footer 'Page N of N' and it has no XObject images
    or annotations.
    """
    try:
        import pypdf
    except ImportError:
        print("⚠️  pypdf not available, skipping blank-page trim")
        return

    reader = pypdf.PdfReader(str(pdf_path))
    writer = pypdf.PdfWriter()
    pages = reader.pages
    n = len(pages)
    keep = n

    for i in range(n - 1, -1, -1):
        page = pages[i]
        text = re.sub(r'\s+', '', page.extract_text() or '')
        footer_only = re.fullmatch(r'Page\d+of\d+', text) is not None
        resources = page.get('/Resources') or {}
        has_images = '/XObject' in resources
        has_annots = bool(page.get('/Annots'))
        if footer_only and not has_images and not has_annots:
            keep = i
        else:
            break

    if keep == n:
        print("   No trailing blank pages")
        return

    for i in range(keep):
        writer.add_page(pages[i])
    with open(pdf_path, 'wb') as f:
        writer.write(f)
    print(f"   Stripped {n - keep} trailing blank page(s): {n} -> {keep}")


def process_pdf_with_metadata(pdf_path: Path):
    """Run the PDF processing script to add metadata and compress."""
    process_script = Path(__file__).parent / 'utils' / 'process_pdf.py'
    
    if not process_script.exists():
        print("⚠️  PDF processing script not found, skipping metadata embedding")
        return
    
    print("🔧 Processing PDF with metadata...")
    try:
        subprocess.run(
            [sys.executable, str(process_script), str(pdf_path), '--quality', 'ebook'],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"⚠️  PDF processing failed: {e}")


async def generate_pdf(quality: str = 'high', wait_time: float = 5.0):
    """Generate PDF from the main TEP paper."""
    
    # Paths — the built static site is the authoritative rendering
    base_dir = Path(__file__).parent.parent
    html_file = base_dir / 'site' / 'dist' / 'index.html'
    
    if not html_file.exists():
        print(f"❌ HTML file not found: {html_file}")
        print("   Run `cd site && npm run build` first.")
        return False
    
    # Load metadata for filename
    metadata = load_citation_metadata()
    version_str = f"v{metadata['version']}-{metadata['codename']}"
    output_name = f"0-TEP-{version_str}.pdf"
    output_pdf = base_dir / output_name
    
    # Select preset based on quality
    # Note: scale controls content zoom (affects page count)
    # device_scale_factor controls pixel density (affects image/text sharpness)
    presets = create_preset_configs()
    if quality == 'maximum':
        # Highest resolution for archival/print quality
        options = presets['high_quality'].copy()
        options['scale'] = 0.72  # Larger content, target <20 pages
        options['device_scale_factor'] = 3.0  # High pixel density for sharpness
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'high':
        options = presets['high_quality'].copy()
        options['scale'] = 0.72
        options['device_scale_factor'] = 2.5
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'print':
        options = presets['print_ready'].copy()
        options['scale'] = 0.72
        options['device_scale_factor'] = 2.0
    else:
        options = presets['web_optimized'].copy()
    
    options['wait_time'] = wait_time
    options['format'] = 'A4'
    options['margin_top'] = '1.2cm'
    options['margin_bottom'] = '1.5cm'  # Increased for footer
    options['margin_left'] = '1cm'
    options['margin_right'] = '1cm'
    
    # Prevent a trailing footer-only blank page: the body's 40px bottom padding
    # overflows the last page's printable area, so remove it at print time and
    # strip trailing margins on the final content elements.
    options['custom_css'] = '''
        @media print {
            body { padding-bottom: 0 !important; }
            main > *:last-child,
            article > *:last-child,
            article section > *:last-child,
            article section div.contact-section { margin-bottom: 0 !important; }
        }
    '''
    
    # Add page numbers footer
    options['display_header_footer'] = True
    options['header_template'] = '<div></div>'  # Empty header
    options['footer_template'] = '''
        <div style="font-size:9px; text-align:center; width:100%; color:#555555; font-family:system-ui,-apple-system,sans-serif; padding-bottom:5mm;">
            Page <span class="pageNumber"></span> of <span class="totalPages"></span>
        </div>
    '''
    
    print(f"\n📄 Generating PDF from: {html_file}")
    print(f"   Quality: {quality}")
    print(f"   Wait time: {wait_time}s (for MathJax rendering)")
    print(f"   Output: {output_pdf}")
    
    async with HTMLToPDFConverter() as converter:
        success = await converter.convert_file(
            str(html_file),
            str(output_pdf),
            options
        )
        
        if not success:
            print("❌ PDF generation failed")
            return False
        
        print(f"✅ PDF generated: {output_pdf}")
        print(f"   Size: {output_pdf.stat().st_size / (1024*1024):.2f} MB")
        print(f"   Version: {version_str}")
        
        # Strip any trailing footer-only blank page before metadata/compression
        strip_trailing_blank_pages(output_pdf)
        
        # Process with metadata
        process_pdf_with_metadata(output_pdf)
        
        print(f"\n✅ Complete! PDF available at: {output_pdf}")
        
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Generate PDF from main TEP theory paper',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate with default settings
  python scripts/generate_site_pdf.py
  
  # Maximum quality (highest resolution)
  python scripts/generate_site_pdf.py --quality maximum --wait-time 10
  
  # High quality with longer wait for MathJax
  python scripts/generate_site_pdf.py --quality high --wait-time 10
        """
    )
    
    parser.add_argument(
        '--quality',
        choices=['maximum', 'high', 'print', 'web'],
        default='high',
        help='PDF quality preset (default: high, maximum for highest resolution)'
    )
    parser.add_argument(
        '--wait-time',
        type=float,
        default=5.0,
        help='Seconds to wait for dynamic content (MathJax) to render'
    )
    
    args = parser.parse_args()
    
    try:
        success = asyncio.run(generate_pdf(
            quality=args.quality,
            wait_time=args.wait_time
        ))
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⚠️  Cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()

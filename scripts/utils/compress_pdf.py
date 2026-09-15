#!/usr/bin/env python3
"""Compress PDF using Ghostscript."""

import subprocess
import shutil
from pathlib import Path


def compress_pdf(input_path: str, output_path: str, quality: str = 'ebook'):
    """Compress a PDF with Ghostscript.

    Returns a dict with:
        original_size: int (bytes)
        compressed_size: int (bytes)
        reduction_percent: float
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    original_size = input_path.stat().st_size

    gs = shutil.which('gs')
    if not gs:
        # Fallback: just copy
        output_path.write_bytes(input_path.read_bytes())
        return {
            'original_size': original_size,
            'compressed_size': original_size,
            'reduction_percent': 0.0,
        }

    quality_settings = {
        'screen': '/screen',
        'ebook': '/ebook',
        'printer': '/printer',
        'prepress': '/prepress',
        'default': '/default',
    }
    pdfsettings = quality_settings.get(quality, '/ebook')

    cmd = [
        gs,
        '-sDEVICE=pdfwrite',
        '-dCompatibilityLevel=1.4',
        f'-dPDFSETTINGS={pdfsettings}',
        '-dNOPAUSE',
        '-dQUIET',
        '-dBATCH',
        f'-sOutputFile={output_path}',
        str(input_path),
    ]

    subprocess.run(cmd, check=True, capture_output=True)

    compressed_size = output_path.stat().st_size
    reduction = (
        (original_size - compressed_size) / original_size * 100
        if original_size > 0 else 0.0
    )

    return {
        'original_size': original_size,
        'compressed_size': compressed_size,
        'reduction_percent': reduction,
    }

import markdown
from playwright.sync_api import sync_playwright
import os
import sys

def convert_md_to_pdf(md_path, pdf_path):
    print(f"Converting {md_path} to {pdf_path}...")
    try:
        with open(md_path, 'r', encoding='utf-8') as f:
            md_content = f.read()

        html_content = markdown.markdown(md_content, extensions=['tables', 'fenced_code', 'toc'])
        
        # Add basic styling for better PDF look
        styled_html = f"""
        <html>
        <head>
        <style>
            body {{ font-family: sans-serif; line-height: 1.6; padding: 40px; }}
            h1, h2, h3 {{ color: #1e3a5f; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            code {{ background-color: #f4f4f4; padding: 2px 4px; border-radius: 4px; }}
            pre {{ background-color: #f4f4f4; padding: 10px; border-radius: 4px; overflow-x: auto; }}
        </style>
        </head>
        <body>
            {html_content}
        </body>
        </html>
        """

        with sync_playwright() as p:
            browser = p.chromium.launch()
            page = browser.new_page()
            page.set_content(styled_html)
            page.pdf(path=pdf_path, format="A4", margin={"top": "1cm", "bottom": "1cm", "left": "1cm", "right": "1cm"})
            browser.close()
        print("Done.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    docs_dir = "docs"
    files_to_convert = ["Optimiser_Manual.md", "SPR_Manual.md", "Optimisation_Workflow_Manual.md"]
    
    for filename in files_to_convert:
        md_file = os.path.join(docs_dir, filename)
        pdf_file = os.path.join(docs_dir, filename.replace(".md", ".pdf"))
        if os.path.exists(md_file):
            convert_md_to_pdf(md_file, pdf_file)
        else:
            print(f"File not found: {md_file}")

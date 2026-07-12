# iolite Optimiser - Workspace Roadmap & TODOs

This file tracks upcoming tasks, design decisions, and roadmap items that span across chat sessions.

---

## 📋 Active Tasks & Decisions

### 1. [ ] Automate PDF Manual Compilation in GitHub Actions
* **Description**: Instead of compiling Markdown documents to PDF locally before a release, configure the GitHub Actions release workflow to compile them dynamically.
* **Proposed Workflow Changes**:
  Add Python, markdown, and Playwright dependencies to `.github/workflows/package-release.yml` to compile PDFs using `convert_docs.py` before creating the release ZIP.
* **Design Decision Needed**:
  * **Approach A (Recommended)**: Delete `.pdf` files from Git, add `docs/*.pdf` to `.gitignore`, and generate them solely on-the-fly inside the GitHub Action for the release `.zip`.
  * **Approach B**: Keep PDFs in Git and have the GitHub Action auto-commit compiled PDFs back to the repository.

---

## 💡 Future Ideas / Backlog
* **[ ]** Track any other minor enhancements or adjustments to the Pulse Train Simulator.

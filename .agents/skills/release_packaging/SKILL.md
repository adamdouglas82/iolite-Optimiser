---
name: Release Packaging and Documentation Compiler
description: Triggered when compiling markdown documentation to PDF, running convert_docs.py, tagging a release, zipping/packaging files, or pushing a new release tag to GitHub.
---

# Release Packaging & Documentation Compilation Workflow

When this skill is triggered, guide the release process systematically following these workflow rules.

---

## 1. Pre-Release Verification & Document Compilation

Before tagging any commit or pushing a release:

* **Audit UI Changes**: Review the Git commits since the last release tag (`git log v<last_version>..HEAD --oneline`) to identify any new features, settings, or UI changes.
* **Verify & Update Documentation**:
  * Check if all identified UI changes are described in their corresponding Markdown manuals in `docs/`.
  * **Spelling Audit**: Check `.md` files for American spelling variants (`optimize`, `analyze`, `normalize`, `color`, `characterize`) and correct to British English (`optimise`, `analyse`, `normalise`, `colour`, `characterise`).
* **Conditional PDF Compilation**:
  > [!IMPORTANT]
  > Execute `python convert_docs.py` ONLY if source `.md` documentation files in `docs/` have been modified. If no Markdown documentation files were changed, do NOT run `convert_docs.py` to prevent unnecessary binary PDF diffs.
  * If `.md` files were updated, run:
    ```powershell
    python convert_docs.py
    ```
  * Stage and include compiled `.pdf` files in the release commit.

---

## 2. Release Commit & Versioning Rules

* **Preserve `VERSION = "dev"`**:
  > [!IMPORTANT]
  > Do **NOT** manually update the `VERSION` constant in `iolite Optimiser.py`. Leave it set to `VERSION = "dev"`. The GitHub Actions release workflow dynamically injects the target version number from the Git tag (`v<version>`) during CI/CD packaging.
* **Draft Multi-Line Release Notes**:
  > [!IMPORTANT]
  > The GitHub Actions workflow (`package-release.yml`) uses `git log -1 --pretty=format:"%B"` to extract the GitHub release notes body. Always commit with a **full multi-line commit message** (subject + section headers + bullet points) so the GitHub release displays complete release notes.
* **Create Release Commit**:
  * Stage all modified codebase files and commit using a file or multi-line string containing the full release notes structure.

---

## 3. GitHub Tagging & CI/CD Automation

* **Create Git Tag**:
  * Tag the release commit locally using semantic versioning:
    ```powershell
    git tag -a v<version> -m "Release v<version>"
    ```
* **Push Release Tag**:
  * Push the tag to GitHub (after explicit user confirmation):
    ```powershell
    git push origin v<version>
    ```
* **CI/CD Workflow Trigger**:
  * Pushing `v<version>` automatically triggers the GitHub Actions `Package Release` workflow, which packages the release zip and publishes the release on GitHub.


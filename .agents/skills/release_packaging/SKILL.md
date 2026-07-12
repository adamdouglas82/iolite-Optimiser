---
name: Release Packaging and Documentation Compiler
description: Triggered when compiling markdown documentation to PDF, running convert_docs.py, tagging a release, zipping/packaging files, or pushing a new release tag to GitHub.
---

# Release Packaging & Documentation Compilation Workflow

When this skill is triggered, you must guide the release process systematically to ensure that compiled PDF manuals are kept in sync with their source markdown files.

---

## 1. Pre-Release Verification & Document Compilation

Before tagging any commit or pushing a release:
1. **Audit UI Changes**: Review the Git commits since the last release tag (using `git log v<last_version>..HEAD --oneline`) to identify any new features, settings, options, or UI layout changes.
2. **Verify & Update Documentation**: 
   * Verify if all identified UI changes and new options are described in their corresponding Markdown manuals in the `docs/` directory.
   * **Automatic Update**: If there are gaps, the agent must **automatically draft the required changes** to the Markdown files (adding sections, updating tables, revising option lists) and present the file edits for user review.
3. **Audit UI Screenshots**: 
   * Check if any modified UI components or panels are represented by screenshots in the manuals (e.g., `![...](docs/images/...)`).
   * Since the agent cannot capture screenshots of the running desktop iolite 4 application, the agent must **explicitly list and alert the user** to any screenshots that need to be re-captured (e.g., `hardware_config_dialog.png` if settings inputs changed) so the user can update the binary image files in Git.
4. **Detect Changes**: Check if any `.md` files in the `docs/` directory have been modified since the last commit or release tag.
5. **Run Compiler**: Execute the Python documentation compiler command to build the latest PDFs:
   ```powershell
   python convert_docs.py
   ```
6. **Commit PDFs**: If running `convert_docs.py` results in new or modified `.pdf` files in the `docs/` folder, stage and commit them:
   * **Commit Message Standard**: `docs: compile manuals to PDF for release <version>`
   * **Important**: Always show the user the proposed commit message and changed files, and get explicit approval before committing.

---

## 2. Release Versioning & Tagging

1. **Verify Version**: Check the target version number (e.g., `1.0.3` or `v1.0.3`).
2. **Update Local Code Version**:
   * Locate `VERSION = "dev"` or the previous version string in [iolite Optimiser.py](file:///c:/Git%20Repositories/iolite%20Optimiser/iolite%20Optimiser.py).
   * Update it to the release version (e.g., `VERSION = "1.0.3"`).
3. **Create Git Tag**:
   * Create the local tag: `git tag -a v<version> -m "Release v<version>"`
   * Always show the user the tag name and message before tagging.

---

## 3. GitHub Push & Release Automation

* **Push Tag**: Ask the user for confirmation to push the tag to GitHub:
  ```powershell
  git push origin v<version>
  ```
* **Release Trigger**: Explain to the user that pushing the tag will automatically trigger the GitHub Actions `Package Release` workflow, which:
  1. Checks out the tagged code (including the newly compiled and committed PDFs).
  2. Updates version metadata.
  3. Archives the code and PDFs into a release `.zip` package.
  4. Publishes it as a GitHub release on the repository.

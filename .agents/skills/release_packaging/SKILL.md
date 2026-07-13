---
name: Release Packaging and Documentation Compiler
description: Triggered when compiling markdown documentation to PDF, running convert_docs.py, tagging a release, zipping/packaging files, or pushing a new release tag to GitHub.
---

# Release Packaging & Documentation Compilation Workflow

When this skill is triggered, you must guide the release process systematically to ensure that compiled PDF manuals are kept in sync with their source markdown files.

---

## 1. Pre-Release Verification & Document Compilation

Before tagging any commit or pushing a release (compilation of PDF manuals by running `convert_docs.py` should only be performed immediately prior to commits, pushes, or releases to avoid redundant compilation of binary PDF files during development):

1. **Audit UI Changes**: Review the Git commits since the last release tag (using `git log v<last_version>..HEAD --oneline`) to identify any new features, settings, options, or UI layout changes.
2. **Verify & Update Documentation**:
   * Verify if all identified UI changes and new options are described in their corresponding Markdown manuals in the `docs/` directory.
   * **Spelling Audit**: Run a grep/search check on the `.md` documentation files for American spelling variants (such as `optimize`, `analyze`, `normalize`, `color`, `characterize`). Correct them to their British English equivalents (`optimise`, `analyse`, `normalise`, `colour`, `characterise`).
   * **Automatic Update**: If there are gaps or spelling inconsistencies, the agent must **automatically draft the required changes** to the Markdown files (adding sections, updating tables, revising option lists) and present the file edits for user review.
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

## 2. Release Commit & Interactive Release Notes

Instead of committing builds separately, consolidate all version updates, manuals, and PDFs into a single Release Commit that triggers the GitHub release notes generator:

1. **Verify Version**: Confirm the target version number (e.g., `1.0.3` or `v1.0.3`).
2. **Update Local Code Version**:
   * Locate the `VERSION` constant in [iolite Optimiser.py](file:///c:/Git%20Repositories/iolite%20Optimiser/iolite%20Optimiser.py) and update it to the target version (e.g., `VERSION = "1.0.3"`).
3. **Draft Release Notes**:
   * Inspect the Git logs since the last tag (`git log v<last_version>..HEAD --oneline`) and summarize all features, bug fixes, and manual changes.
   * Write a structured, multi-line release description draft into a temporary file in the repository root: `release.md`.
4. **User Refinement**:
   * Alert the user that the draft has been written to `release.md` and invite them to open, edit, and format it in their editor.
   * **Wait** for the user's confirmation that they are finished editing.
5. **Create Release Commit**:
   * Read the final contents of `release.md`.
   * Stage all release changes (updated code, markdown manuals, and compiled PDFs).
   * Commit the changes using the contents of `release.md` as the commit message.
   * Delete the temporary `release.md` file from the repository.
6. **Apply Git Tag**:
   * Tag the release commit locally: `git tag -a v<version> -m "Release v<version>"`
   * Show the user the tag details before running the command.

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

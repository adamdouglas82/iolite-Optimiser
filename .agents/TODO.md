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

### 2. [ ] Unload Auto-Imported Laser Log Files in `import_and_unload_data`

* **Description**: Certain mass spectrometer file formats automatically import an accompanying laser log file into iolite when `data.importData(norm_file_path)` is called. Currently, `import_and_unload_data` only unloads the mass spec file via `data.unloadFile(...)`.
* **Current API Status**:
  * `data.importedLogs()` returns a list of `FileMetadataPyInterface` objects (with `.filePath()` and `.fileName()`) for all loaded laser log files.
  * In `iolite_core.dll`, `DataManager` has `DataManager::unloadLaserLogFile(QString filePath)`, but `IolitePythonInterface` (`data`) currently only exposes `data.unloadFile(filePath)` (which calls `DataManager::unloadFileData(filePath)` for mass spec files only).
  * A feature request has been submitted to the iolite software team to expose laser log unloading in the `data` Python API.
* **Implementation Options (When Resumed)**:
  * **Option A (Official Python API)**: Once added by the iolite team (e.g. `data.unloadLog(path)` or updated `data.unloadFile(path)`), record `logs_before = data.importedLogs()` prior to `data.importData(...)` and unload any newly added laser logs in `logs_after` right after `data.unloadFile(...)`.
  * **Option B (UI Fallback if needed earlier)**: Locate iolite's `FilesBrowser` widget (`objectName == "FilesBrowser"`) and its `logsView` table (`objectName == "logsView"`), select the newly imported laser log row(s), ensure `filesView` does not have focus, and invoke `FilesBrowser.unload()` (which calls `DataManager::unloadLaserLogFile(filePath)` directly without a confirmation dialog), with `data.resetSession()` as a fallback when `not files_before and not logs_before`.

---

## 💡 Future Ideas / Backlog
* **[ ]** Track any other minor enhancements or adjustments to the Pulse Train Simulator.

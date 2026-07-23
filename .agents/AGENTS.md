# iolite Optimiser - Workspace Agent Guidelines

These rules and guidelines help Antigravity and other AI coding assistants understand the constraints, architecture, and workflow requirements for the **iolite Optimiser** project.

---

## 1. Environment & Threading Constraints (iolite / PythonQt)

* **Qt Integration**: The optimiser runs as a plugin inside **iolite 4**, using an embedded Python environment integrated with Qt via **PythonQt**.
* **Thread Safety & QTimer Restriction**:
  > [!IMPORTANT]
  > PythonQt / iolite does not allow `QTimer.singleShot` to be called from a non-Qt thread.
  > If you need to perform background/asynchronous tasks (e.g., network calls, checking for updates), you must:
  > 1. Start the task in a background Python thread (`threading.Thread`).
  > 2. Write the result to a thread-safe shared list or container (e.g., `result_holder`).
  > 3. Use a main-thread polling `QTimer` (running every 100 ms or similar) to check for updates and process the result, stopping itself when done.
* **Qt Property & Method Access**:
  > [!NOTE]
  > In PythonQt, Qt widget properties/methods (like `.value()`, `.currentText`, `.count`) can dynamically resolve as either callable methods or direct attribute properties depending on the environment version.
  > * **Always use the safe retrieval helper**: Call `self._get(widget, 'property_or_method_name')` (e.g., `self._get(self.spin_spr_prom, 'value')` instead of `self.spin_spr_prom.value()`).
  > * **Always use the safe signal blocking helper**: Call `self._block_signals(widget, True/False)` instead of `widget.blockSignals(True/False)`.
* **Explicit Base Class Delegation (No `super()`)**:
  > [!WARNING]
  > Standard Python `super()` calls inside overridden Qt widget event methods (like `keyPressEvent`) can fail or misbehave in PythonQt.
  > * **Explicitly invoke the base class methods**: e.g., use `QTableWidget.keyPressEvent(self, event)` instead of `super().keyPressEvent(event)`.
* **Platform/Environment URL Opening**:
  > [!CAUTION]
  > Do **NOT** use Python's standard `webbrowser.open()` to launch URLs. It causes access violation crashes (`0xC0000005`) in iolite's embedded Python 3.13 via `PyOS_AfterFork`/`ShellExecuteW` when called within the Qt event loop.
  > * **Always use native Qt URL handling**: Use `QDesktopServices.openUrl(QUrl(url))`.

---

## 2. Gitflow & Git Operations

* **Branching Strategy (Gitflow)**:
  * This repository follows a Gitflow/feature-branch workflow.
  * Do not make or push commits directly to the `main` branch unless explicitly instructed.
  * Always check current branch status (`git status`, `git branch`) before editing or proposing commits.
* **Review & Confirmations**:
  * **Before Pushing**: Always check with the user and get explicit approval before pushing any changes to GitHub.
  * **Commit Messages**: Always show the user the proposed commit message before committing.
* **Release Versioning & Packaging Workflow**:
  * **Keep `VERSION = "dev"`**: Do not manually change `VERSION` in `iolite Optimiser.py`. The GitHub Actions CI/CD workflow dynamically handles versioning based on the Git release tag (`v<version>`).
  * **Release Commit Title Format**: Always format the top-line title of release commits as `Release v<version>: <summary>` (e.g. `Release v1.0.6: pulse train simulator updates and code cleanup`). Do not prefix with `release: release...`.
  * **Conditional PDF Compilation**: Only run `convert_docs.py` if source `.md` files in `docs/` have been modified. If no Markdown documentation files were changed, do not run `convert_docs.py` or commit unneeded PDF binary diffs.
  * **Release Tagging**: Always tag release commits with semantic versioning (`git tag -a v<version> -m "Release v<version>"`) and push the tag (`git push origin v<version>`) after explicit user approval to trigger GitHub Actions release packaging.

---

## 3. Session Continuity & Memory

* **Task Tracking**:
  * At the start of a session, check if there is an active `task.md` or task list.
  * If a task is ongoing, read the recent history or status files to pick up where the previous session left off.
  * Update `task.md` or write status summaries regularly to assist the next session.

---

## 4. Coding Style & Language Standards

* **British English Localisation**:
  * All user-facing UI labels, documentation manuals, and code implementations (variable names, method names, comments) must use **British English** spelling.
  * **Key words**: Use `optimisation` / `optimise`, `analyse` / `analysis`, `normalise` / `normalisation`, `colour` (e.g. `run_optimisation`, `analyse_washout_peaks`, etc.). Do not use `optimization`, `optimize`, `analyze`, or `color`.
* **Plain Text Mathematical Formatting**:
  * Do **NOT** use LaTeX math delimiters (such as `\(...\)`, `\[...\]`, or `$$...$$`) in conversation responses, commit messages, or markdown documents. The user's chat window does not render LaTeX correctly.
  * Always format mathematical expressions, units, and formulas in clear, plain text (e.g. `11.08 µm`, `325 Hz`, `299.00 µm/s`, `Speed = RepRate * (SpotSize - Overlap)`).
* **Markdown Formatting & Linting**:
  * **List Consistency**: Always use asterisk (`*`) for all unordered list items (including nested items) in markdown files (such as `session_history.md`, `walkthrough.md`, `implementation_plan.md`, and `task.md`) to prevent list style mixed warnings.
  * **Blank Lines**: Always surround headings and lists with exactly one blank line (preventing `blanks-around-headings` and `blanks-around-lists` warnings).
  * **File Headings**: Always start markdown files with a top-level `#` H1 heading as the first line of the file (preventing `first-line-h1` warnings).

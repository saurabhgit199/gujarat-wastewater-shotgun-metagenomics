
C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>git commit -m "Add comprehensive README with full project explanation and findings"
On branch main
Your branch is ahead of 'origin/main' by 2 commits.
  (use "git push" to publish your local commits)

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
        modified:   .gitignore

Untracked files:
  (use "git add <file>..." to include in what will be committed)
        reports/
        requirements.txt
        scripts/.Rhistory

no changes added to commit (use "git add" and/or "git commit -a")

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>git push
Enumerating objects: 27, done.
Counting objects: 100% (27/27), done.
Delta compression using up to 12 threads
Compressing objects: 100% (24/24), done.
error: RPC failed; HTTP 408 curl 22 The requested URL returned error: 408
send-pack: unexpected disconnect while reading sideband packet
Writing objects: 100% (25/25), 36.20 MiB | 290.00 KiB/s, done.
Total 25 (delta 5), reused 0 (delta 0), pack-reused 0 (from 0)
fatal: the remote end hung up unexpectedly
Everything up-to-date

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>powershell -Command "Invoke-WebRequest -Uri 'https://raw.githubusercontent.com/saurabhgit199/gujarat-wastewater-shotgun-metagenomics/main/README.md' -OutFile 'README_current.md'"

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>cd "C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github"

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>notepad README.md

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>git add README.md

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>git commit -m "Add comprehensive README with full project explanation"
On branch main
Your branch is ahead of 'origin/main' by 2 commits.
  (use "git push" to publish your local commits)

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
        modified:   .gitignore

Untracked files:
  (use "git add <file>..." to include in what will be committed)
        README_current.md
        reports/
        requirements.txt
        scripts/.Rhistory

no changes added to commit (use "git add" and/or "git commit -a")

C:\Users\SAURABH\OneDrive\Desktop\lotus falk work\metagene-github>git push
Enumerating objects: 27, done.
Counting objects: 100% (27/27), done.
Delta compression using up to 12 threads
Compressing objects: 100% (24/24), 

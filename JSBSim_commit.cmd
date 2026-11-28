@rem Before a Visual Studio build of the library: src\JSBSim_commit.h names the commit it is built from, rewritten only when that changes.
set JSBSIM_COMMIT=
if exist .git for /f %%c in ('git rev-parse --verify HEAD') do set JSBSIM_COMMIT=%%c
if not defined JSBSIM_COMMIT (echo JSBSim names the commit it is built from, and git rev-parse --verify HEAD names none here & exit /b 1)
findstr /c:"%JSBSIM_COMMIT%" src\JSBSim_commit.h >nul 2>&1 || echo #define JSBSIM_COMMIT "%JSBSIM_COMMIT%">src\JSBSim_commit.h
findstr /c:"%JSBSIM_COMMIT%" src\JSBSim_commit.h >nul

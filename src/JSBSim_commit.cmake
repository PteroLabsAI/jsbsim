# At configure and, in a git checkout, at every build of libJSBSim
# (src/CMakeLists.txt): HEADER names the commit of the git checkout SOURCE,
# rewritten only when that changes, or stays as the source package it is in
# (PKG-INFO in SOURCE) was made with.
cmake_minimum_required(VERSION 3.18)
if(EXISTS ${SOURCE}/.git)
  execute_process(COMMAND git rev-parse --verify HEAD
                  WORKING_DIRECTORY ${SOURCE}
                  OUTPUT_VARIABLE JSBSIM_COMMIT
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE JSBSIM_COMMIT_FAILED)
  if(JSBSIM_COMMIT_FAILED)
    message(FATAL_ERROR "git rev-parse --verify HEAD failed in ${SOURCE}: JSBSim names the commit it is built from")
  endif()
  file(CONFIGURE OUTPUT ${HEADER}
       CONTENT "#define JSBSIM_COMMIT \"@JSBSIM_COMMIT@\"\n" @ONLY)
elseif(NOT (EXISTS ${SOURCE}/PKG-INFO AND EXISTS ${HEADER}))
  message(FATAL_ERROR "${SOURCE} is neither a git checkout nor a source package holding ${HEADER}: JSBSim names the commit it is built from")
endif()

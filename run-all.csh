#! /bin/tcsh -f
# An entry point for running all of the carma tests. It creates a run directory
# copies in the executables and then runs all of the test executables.
#
#
# Usage
#  run-carma.csh

# Environment Variables
#   CARMA_BUILD [carma]
#     The subdirectory in which the build was performed.

# By default, run the single column model
if (! $?CARMA_BUILD ) then
  setenv CARMA_BUILD carma
endif

if (! $?CARMA_CASE ) then
  setenv CARMA_CASE $CARMA_BUILD
endif

if (! $?CARMA_THREADS ) then
  setenv CARMA_THREADS 1
endif

#set runtgt=CARMATEST.exe

#if ($# == 1) then
#  set runtgt="$1"
#endif

set blddir=build/$CARMA_BUILD
set rundir=run/$CARMA_CASE
set testdir=tests

# Create a directory for the build.
mkdir -p $rundir

# Copy the executable to the run directory.
cp $blddir/*.exe $rundir


# Prepare for multiple threads, assuming Intel Compiler.
setenv OMP_NUM_THREADS $CARMA_THREADS
setenv KMP_STACKSIZE 128M

# Execute the tests.
cd $rundir

echo `ls -1 *.exe`
foreach runtgt (`ls -1 *.exe`)
  echo "  ** Started $runtgt at `date` **"
  ./$runtgt || echo '  *** Run Failed ***' && exit -1
  echo "  ** Finished at `date` **"
  echo ""
end



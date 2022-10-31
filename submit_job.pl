#!/usr/bin/perl -w

use Cwd;

$package = "pp13TeV_gg2qq";
$maindir = getcwd();

$groupnum = 0;
$progname = "mainEx04";

$rundir = "${maindir}/run_${package}_grp${groupnum}";
mkdir $rundir;

$jobname = sprintf("${package}_grp%03d",${groupnum});

for ($irun=0; $irun<10; $irun++){

#	sleep 1;

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
	print FILE "Getenv = true\n";
	print FILE "Priority = +30\n";
	print FILE "Executable = jobscript\n";
	print FILE "JobBatchName = ${jobname}\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
	print FILE "Queue\n";
	close(FILE);

	$seed = $irun+21;

	open(FILE, ">${progname}.cfg");
	print FILE "Main:numberOfEvents = 10000000\n";
	print FILE "Main:timesAllowErrors = 3\n";
	print FILE "Random:setSeed = on\n";
	print FILE "Random:seed = ${seed}\n";
	print FILE "Init:showChangedSettings = on\n";
	print FILE "Init:showChangedParticleData = on\n";
	print FILE "Next:numberCount = 1000\n";
	print FILE "Next:numberShowInfo = 1\n";
	print FILE "Next:numberShowProcess = 1\n";
	print FILE "Next:numberShowEvent = 1\n";
	print FILE "Beams:idA = 2212\n";
	print FILE "Beams:idB = 2212\n";
	print FILE "Beams:eCM = 13000.\n";
	print FILE "HardQCD:gg2qqbar = on\n";
	print FILE "PartonLevel:MPI = off\n";
    print FILE "PartonLevel:ISR = off\n";
	print FILE "PartonLevel:FSR = off\n";
	print FILE "PhaseSpace:pTHatMin = 5\n";

	close(FILE);

#	$seednum = int(rand(100000000));

	$grpdir = sprintf("%s/out_%s_grp%03d",$maindir,$package,$groupnum);

#	$randnum = int(rand(60));

	open(FILE, ">jobscript");
	print FILE "#!/bin/bash\n";
	print FILE "mkdir -p $grpdir\n\n";

	print FILE "cp -av ${maindir}/${progname} .\n";

	print FILE "./${progname} 5\n";

	$outputfile = sprintf("%s/outfile_%s_grp%03d_%05d.root",$grpdir,$package,$groupnum,$irun);
	print FILE "mv Pythia8_event.root $outputfile\n";

	print FILE "rm -rf ${progname}\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";
}


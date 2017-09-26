#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/lib";
use AnaMethod;
use Data::Dumper;

my ($input,$outdir,$plot,$env,$move, $monitorOption,$java, $help, $artifactModes, $qsubMemory, $pymonitor, $python, $Rscript);

GetOptions(
	"input:s" => \$input,
	"env:s"=> \$env,
	"outdir:s" => \$outdir,
	"move:s" => \$move,
	"m:s" => \$monitorOption,
	"help|?" => \$help,
	"qsubMemory:s" => \$qsubMemory,
	"pymonitor:s" => \$pymonitor,
	"python:s" => \$python,
	"Rscript:s" => \$Rscript,
	"plot" => \$plot,
);
my $usage = <<USE;
Usage:
description: SV Signature analysis
author: LuoShizhi, luoshizhi\@genomics.cn
date: 2017-08-17
usage: perl $0 [options]
	Common options:
	-input*		<str>	breakpoint file dir. keyword\tbreakpoint_file_dir\tdependence.sh(optioanl)
	-env <str>   export environment variable[\$Bin/environment.sh]
				 PATH/MANPATH/LD_LIBRARY_PATH/CFLAGS/LDFLAGS/C_INCLUDE_PATH/CPLUS_INCLUDE_PATH/LIBRARY_PATH/CPATH/R_LIBS
	-outdir*		<str>	outdir.[./]
	-move		<str>	if this parameter is set,final result will be moved to it from output dir.
	-m		<str>	monitor options. will create monitor shell while defined this option
	-qsubMemory	<str>	Job memory.format:step1_mem,step2_mem[1G,8G,8G]
	-help|?			print help information
	-plot  plot a pic

	Database options:
	-fa		<tr>	genome fa.[\$Bin/database/hg19.fa]
	Software options:
	-pymonitor	<str>	monitor path [\$Bin/bin/monitor]
	-Rscript	<str>	Rscript path [\$Bin/bin/Rscript]

e.g.:
	perl $0 -input breakpoint.list -outdir ./outdir
USE

die $usage unless ($input && $outdir);
$qsubMemory ||= "1G,4G,8G";
my @qsubMemory = split /,/,$qsubMemory;
$qsubMemory[0] ||= "1G";
$qsubMemory[1] ||= "4G";
$qsubMemory[1] ||= "8G";
$outdir ||= "./";
mkpath $outdir;
$outdir = File::Spec->rel2abs($outdir);

if($move){
	$move = File::Spec->rel2abs($move);
	mkpath $move;
}
$Rscript ||= "$Bin/bin/Rscript";
$env ||="$Bin/environment.sh";
$monitorOption ||="-P common -q bc.q -p test";

my ($breakpoint,$dependent) = &ReadInfo($input);
my %breakpoint = %$breakpoint;
my %dependent = %$dependent;


my ($shell, $process, $list)=("$outdir/shell/", "$outdir/process/", "$outdir/list");
mkpath $shell;mkpath $process;mkpath $list;
my $dependence = "$list/dependence.txt";

open TXT, ">$dependence" or die $!;
my $out  = ($move)?$move:$process;

#print Dumper %breakpoint;
#die;
foreach my $breakfile (keys %breakpoint){
	my $process_t = "$process/$breakfile"; mkpath $process_t;
	my $shell_t = "$shell/$breakfile" ; mkpath $shell_t;
	###step1 sv_clust
	my $sv_clust = "$shell_t/step1.sv_clust.sh";
	my $content ="source $env&&\\\n";
	$content .="cd $process_t  &&\\\n";
	$content .="$Rscript $Bin/script/sv_clust.r $breakpoint{$breakfile}";
	AnaMethod::generateShell($sv_clust,$content);
	if (exists $dependent{$breakfile}) {
		print TXT "$dependent{$breakfile}\t$sv_clust:$qsubMemory[0]\n";
	}
	
	###step2 run_Signatures
	my $run_Signatures="$shell_t/step2.run_Signatures.sh";
	$content ="source $env&&\\\n";
	$content .="cd $process_t  &&\\\n";
	$content .="$Bin/script/run_Signature.sh &&\\\n";
	$content .="$Rscript $Bin/script/plot.nmf.R $process_t/output/parameter_index.txt $process_t/rank_vs_nmf_metrics.pdf";
	AnaMethod::generateShell($run_Signatures,$content);
	print TXT "$sv_clust:$qsubMemory[0]\t$run_Signatures:$qsubMemory[0]\n";
	###step3 

	if (defined $plot) {
		my @content;
		my @files=glob("$process_t/output/*processes.txt");
		my $plotsh="$shell_t/step3.sh";
		#foreach my $file (@files){	
			#push @content,"$Rscript $Bin/script/plot_signature.r $file sv_signature.$file";
		#}
		#$content=join"&&\\\n",@content;
		$content="for file in $process_t/output/*processes.txt\n";
		$content .="do\n";
		$content .="$Rscript $Bin/script/plot_signature.r \$file \$file.sv_signature\n";
		$content .="done";
		AnaMethod::generateShell($plotsh,$content);
		print TXT "$run_Signatures:$qsubMemory[0]\t$plotsh:$qsubMemory[0]\n";
	}
}
if(defined $pymonitor && defined $monitorOption){
	`echo "$pymonitor $monitorOption -i $dependence" >$list/qsub.sh`;
}

close TXT;


sub ReadSampleInfo {
        my ($file) = @_;
        my (%hashSample,%hashDepend,%hashlength);
        open IN, "$file" or die $!;
        while (<IN>) {
                chomp;
                next if(/^\s*$/);
                s/\s*$//;
                s/^\s*//;
                my @tmp = split /\t+/;
                $hashSample{$tmp[0]}=$tmp[1];
				$hashlength{$tmp[0]}=$tmp[2];
                $hashDepend{$tmp[0]}=$tmp[3] if(@tmp >= 3);
        }
        close IN;
        return (\%hashSample,\%hashDepend,\%hashlength);
}
sub ReadInfo {
        my ($file) = @_;
        my (%hashbreakpoint,%hashDepend);
        open IN, "$file" or die $!;
        while (<IN>) {
                chomp;
                next if(/^\s*$/);
                s/\s*$//;
                s/^\s*//;
                my @tmp = split;
                $hashbreakpoint{$tmp[0]}=$tmp[1];
                $hashDepend{$tmp[0]}=$tmp[2] if(@tmp >= 3);
        }
        close IN;
        return (\%hashbreakpoint,\%hashDepend);
}

sub Readpair2 {
        my ($file) = @_;
        my ($control,$treatment,%pair,%C,%T,%T_N);
        open IN, "$file" or die $!;
        while (<IN>) {
                next if(/^\s*$/);
                chomp;
                s/\s*$//;
                s/^\s*//;
                next if /^\s+#/;
                if(/(\S+)\t(\S+)/){
                        $control = $1;
                        $treatment = $2;
                }
                $T_N{$treatment}=$control;
        }
        close IN;
        return (%T_N);
}


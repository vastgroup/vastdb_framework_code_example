#!/usr/bin/perl
# Original name: get_tisAS_v6.pl
# Slightly modified from the version in https://github.com/vastdb-pastdb/pastdb
# group file: Sample_name\tGroup\tGr_Excluded1,Gr_Excluded2
use strict;
no strict "refs";
use warnings;
use Getopt::Long;

#### General variables
my $min_dPSI=15;
my $min_dPSI_glob=25;
my $min_dPSI_strict=5;
my $min_range=0; # to define AS among the sample with coverage
my $print;
my $test_type="EX";
my $Nt=10;
my $noVLOW;
my $p_IR;
my $median;
my $test_tis;
my $min_rep=1;
my $Q="O[WK]\,.+?\,.+?\,.+?\,.+?\@";
my $strict;
my $groups;
my $use_truly;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_dPSI=i" => \$min_dPSI,
			  "min_dPSI_glob=i" => \$min_dPSI_glob,
			  "min_dPSI_strict=i" => \$min_dPSI_strict,
                          "min_range=i" => \$min_range,
			  "event_type=s" => \$test_type,
			  "N_groups=i" => \$Nt,
			  "N=i" => \$Nt,
                          "groups=s" => \$groups,
                          "g=s" => \$groups,
			  "use_truly=s" => \$use_truly,
			  "print" => \$print,
			  "strict" => \$strict,
			  "median" => \$median,
			  "test_tis=s" => \$test_tis,
			  "min_rep=i" => \$min_rep,
                          "help" => \$helpFlag,
			  "p_IR" => \$p_IR,
                          "noVLOW" => \$noVLOW
    );

if (!defined($ARGV[0]) || $helpFlag){
    die "\nUsage: Get_Tissue_Specific_AS.pl INCLUSION_LEVELS_FULL-root.tab -g Tissue_groups [options]

Identifies tissue-specific events by assessing non-overlapping PSI of groups

[General options] 
        --min_dPSI i             Minimum delta PSI of the average against all other averages (default 15)
        --min_dPSI_glob i        Minimum delta PSI of the average against the average of all other samples (default 25)
        --min_range i            Minimum PSI range for the event to be considered AS among groups (default 0)
        --groups/-g file         Groups file (if not provided, each sample is a group)
        --use_truly i            Uses a truly file to define AS, instead of range (default OFF)
        --N_groups/-N i          Minumum number of groups with coverage, *including* target tissue (default 10)
        --test_tis               Comma-separated list of issues for which it will make a table with dPSIs (default OFF)
        --event_type s           Type of AS event to be tested (EX, IR, Alt3, Alt5, BS; default EX)
        --noVLOW                 Does not use samples with VLOW coverage (default OFF)
        --p_IR                   Applies the binomial test filter for IR (default OFF)
        --strict                 Requires all sub-samples to also be non-overlapping within PSI min_dPSI_strict (default OFF)
        --min_dPSI_strict i      Minimum dPSI between each of the subsamples.
        --median                 Uses PSI median of group, instead of average (default OFF)
        --min_rep                Minimum number of replicates with coverage for a group to be considered (default 1)
        --print                  Print an output file with the tissue-specific events (default OFF)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

my $v;
($v)=$ARGV[0]=~/INCLUSION_LEVELS[\-\_].+?\-(.{3}.+?)[\-\.]/; # to allow subsets of V
if (!$v){
   ($v)=$ARGV[0]=~/\-(.+?)\./;
}

my %truly;
if ($use_truly){
    open (TRULY, $use_truly) || die "Truly";
    <TRULY>;
    while (<TRULY>){
	chomp;
	my @t=split(/\t/);
	$truly{$t[1]}=1 if $t[9] eq "SuperAS" || $t[9] eq "TrulyAS";
    }
    close TRULY;
}


my $output_file="TS_AS-$v-N$Nt-dPSI$min_dPSI-dPSI_glob$min_dPSI_glob-Range$min_range-Rep$min_rep-$test_type";
$output_file.="-noVLOW" if $noVLOW;
$output_file.="-p_IR" if $p_IR;
$output_file.="-Median" if $median;
$output_file.="-Strict" if $strict;
$output_file.="$min_dPSI_strict" if $strict;

if (defined $print){
    open (O, ">$output_file.tab");
    print O "GENE\tEventID\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tTissue\tDirection\tN_Rep\tTotal_N_others\tAv_Tis\tAv_Others\tdPSI\t".
       "SD_Tis\tSD_Others\tMin_Others\tMax_Others\tValid_Tis\n";
}

### Creates the files for each test tissues
my @TEST_TIS;
if (defined $test_tis){
    @TEST_TIS = split(/\,/,$test_tis);
}
my %test_tis_OK;
foreach my $tis (@TEST_TIS){
    my  $handler="dPSIs_$tis";
    open ($handler, ">dPSIs_$tis-$v-N$Nt-Rep$min_rep-$test_type.tab");
    print $handler "GENE\tEventID\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tN_Rep\tTotal_Others\tAv_$tis\tAv_Others\tdPSI_$tis-Others\t".
	"SD_Tis\tSD_Others\tMin_Others\tMax_Others\tValid_samples\n";
    $test_tis_OK{$tis}=1;
}


my %group;
my %temp_excl;
my %excl;
if ($groups){
    open (GROUPS, $groups) || die "Can't open group file";
    while (<GROUPS>){
	s/\r//g;
	s/\"//g;
	s/ //g;
	chomp;
	my @t=split(/\t/);
	$group{$t[0]}=$t[1] if $t[1];
	$group{$t[0]}=$t[0] if !$t[1];
	
        if (defined $t[2]){
            if (defined $group{$t[0]}){
            if (defined $temp_excl{$group{$t[0]}}){
	    die "Inconsistent exclusion set for the same group ($group{$t[0]}: $temp_excl{$group{$t[0]}})\n" if $temp_excl{$group{$t[0]}} ne $t[2] && $temp_excl{$group{$t[0]}} && $group{$t[0]} ne "EXCLUDE";
        }}}
	$temp_excl{$group{$t[0]}}=$t[2]; # comma separated (per group), to check
	my @excluded;
        if (defined $t[2]){
            @excluded = split(/\,/,$t[2]);  
        }
	foreach my $temp (@excluded){
	    $excl{$group{$t[0]}}{$temp}=1; # real test of the actual excluded tissue
	}
    }
    close GROUPS;
}

if ($ARGV[0]=~/\.gz$/){
   open (I, "gunzip -c $ARGV[0] | ") || die "Cannot open INCL $ARGV[0]";
} else {
   open (I, $ARGV[0]) || die "Cannot open INCL $ARGV[0]"; # INCL table
}
my $head=<I>;
chomp($head);
my @head=split(/\t/,$head);
my %val; my %Q;
for my $i (6..$#head){
    if ($i%2==0){
	$val{$head[$i]}=$i;
	$Q{$head[$i]}=$i+1;
	print "  No group identified for $head[$i]\n" if !$group{$head[$i]} && $groups;
	$group{$head[$i]}=$head[$i] if !$groups;
    }
}

my %tallyTIS_AS;
my %tally; # keeps the count of tissue-specific  
my (%for_matV, %KD, %TOTOKs, %SUMOKs, %sum, %tot, %max, %min, %OK, %vals);    
my (%med, %av, %SD, %TOTAL, %min_ev, %max_ev);
while (<I>){
    s/VLOW/N/g if $noVLOW;
    chomp;
    my @t=split(/\t/);
    
    my $ev=$t[1];
    my $actualtis=0;
    my $pre_line=join("\t",@t[0..5]);
    my $valid_samples="";

    my $type="";
    $type="IR" if $t[5]=~/IR/;
    $type=$t[5] if $t[5]=~/Alt/;
    $type="EX" if $t[1]=~/EX/;
    $type="BS" if $t[5]=~/BS/; # future circRNAs
    
#    next if ($t[3]!=0 && $t[5]=~/Alt/); # only the shorter version (i.e. 1 per Alt event)
    next if $t[3]==0;
    next if $type!~/$test_type/i;
    next if !$truly{$ev} && $use_truly; # only proper AS (but then check below)

    for my $i (6..$#t){
	if ($i%2==0 && $t[$i+1]=~/$Q/){ # it analyzes per column, rather than before
	    my $tis=$head[$i];
	    my $kill_pIR = "";
	    if ($t[5]=~/IR/ && $p_IR){
                my ($temp_p)=$t[$Q{$tis}]=~/\,.+?\,.+?\,.+?\,(.+?)\@/; # previously tis not i
                $kill_pIR = 1 if $temp_p < 0.05;
            }
	    my $AT;
            if ($group{$tis}){$AT=$group{$tis};}
	    next if !$AT;
	    next if $AT eq "EXCLUDE"; # if the group is EXCLUDE, exclude.
	    next if $kill_pIR;
	    parse_sample($ev,$tis,$AT,@t);
	}
    }
####
    
#DO the total tissue (group, in fact)
    foreach my $tis (sort (keys %{$sum{$ev}})){
	if ($tot{$ev}{$tis} >= $min_rep){ # i.e. the number of samples for the group is >= min_rep
	    $av{$ev}{$tis}=sprintf("%.2f",$sum{$ev}{$tis}/$tot{$ev}{$tis});
	    $med{$ev}{$tis}=median($ev,$tis);
	    $SD{$ev}{$tis}=std_dev($ev,$tis,$av{$ev}{$tis}) if $tot{$ev}{$tis}>1; # back in V5
	    $SD{$ev}{$tis}="NA" if $tot{$ev}{$tis}==1;
	    $av{$ev}{$tis}=$med{$ev}{$tis} if $median; ## to compare medians
            if (defined $min_ev{$ev}){
	         $min_ev{$ev}=$av{$ev}{$tis} if $av{$ev}{$tis} <= $min_ev{$ev};
            } else {
                 $min_ev{$ev}=$av{$ev}{$tis};
            }
            if (defined $max_ev{$ev}){
	         $max_ev{$ev}=$av{$ev}{$tis} if $av{$ev}{$tis} >= $max_ev{$ev};
            } else {
                 $max_ev{$ev}=$av{$ev}{$tis};
            }  
	    $TOTAL{$ev}++;
	    $valid_samples.="$tis,";
	}
    }
    chop($valid_samples);
    
    if (defined $TOTAL{$ev}){
    if ($TOTAL{$ev}>=$Nt){ # Nt is the total number of groups (including the one to be compared)
	next if $max_ev{$ev}-$min_ev{$ev} < $min_range; # AS within the group.
	
	foreach my $tis1 (sort (keys %{$av{$ev}})){
	    my ($real_total, $pos, $neg)=(0,0,0);
	    my $sum_others=0;
	    my @vals_others=();

	    foreach my  $tis2 (sort (keys %{$av{$ev}})){
		next if $excl{$tis1}{$tis2}; # 27/03/18

		if ($tis1 ne $tis2){
		    $real_total++; # the total number of tissues matched against tis1
		    $sum_others+=$av{$ev}{$tis2};
		    push(@vals_others,$av{$ev}{$tis2});
		    if ($strict){ # no overlap for any sub-sample
			$pos++ if ($av{$ev}{$tis1}-$av{$ev}{$tis2} >= $min_dPSI && $min{$ev}{$tis1}-$max{$ev}{$tis2}>=$min_dPSI_strict);
			$neg++ if ($av{$ev}{$tis1}-$av{$ev}{$tis2} <= -1*$min_dPSI && $max{$ev}{$tis1}-$min{$ev}{$tis2} <= -$min_dPSI_strict);
		    }
		    else {
			$pos++ if ($av{$ev}{$tis1}-$av{$ev}{$tis2} >= $min_dPSI);
			$neg++ if ($av{$ev}{$tis1}-$av{$ev}{$tis2} <= -1*$min_dPSI);
		    }
		}
	    }
	    
	    if ($real_total >= $Nt-1){ #total minus tis1 but after filtering for excluded
		$tallyTIS_AS{$tis1}++; # keeps the count of valid tests
		my $dir="";
		$dir="DOWN" if $neg == $real_total; 
		$dir="UP" if $pos == $real_total;
		my $av_others=sprintf("%.2f",$sum_others/$real_total); # the average of the compared tissues
		my $min_others=(sort{$a<=>$b}(@vals_others))[0];
		my $max_others=(sort{$b<=>$a}(@vals_others))[0];
		my $dPSI_tis1=sprintf("%.2f",$av{$ev}{$tis1}-$av_others);
		my $handler="dPSIs_$tis1";
		my $SD_others=std_dev2(@vals_others); ### Calculates the SD for "others" (new in V5)
		
		print $handler "$pre_line\t$tot{$ev}{$tis1}\t$real_total\t$av{$ev}{$tis1}\t$av_others\t$dPSI_tis1\t".
		    "$SD{$ev}{$tis1}\t$SD_others\t$min_others\t$max_others\t$valid_samples\n" if $test_tis_OK{$tis1};
		
		if ($dir && abs($dPSI_tis1)>=$min_dPSI_glob){
		    $tally{$tis1}{$dir}++; # keeps the count of tissue-specific
                    if (defined $print){
		         print O "$pre_line\t$tis1\t$dir\t$tot{$ev}{$tis1}\t$real_total\t$av{$ev}{$tis1}\t$av_others\t$dPSI_tis1\t".
			   "$SD{$ev}{$tis1}\t$SD_others\t$min_others\t$max_others\t$valid_samples\n";
                    }  
		}
	    }
	}
    }
    } # if defined $TOTAL{$ev}
}

print "$v: min_dPSI>$min_dPSI, min_dPSI_glob>$min_dPSI_glob, N_groups $Nt, min_rep $min_rep";
print ", noVLOW" if $noVLOW;
print ", TrulyAS" if $use_truly;
print ", Min_range $min_range" if $min_range;
print ", Strict ($min_dPSI_strict)" if $strict;
print ", Median" if $median;
print "\n";

print "SAMPLE\tTOTAL\tUP\tDOWN\tTot_EV\tPerc_TS-EV\tPerc_UP\tPerc_DOWN\n"; 
my %data_print;
my @order;
foreach my $tis (sort (keys %tallyTIS_AS)){ # which eq tallyTIS_ALL
    $tally{$tis}{UP}=0 if !$tally{$tis}{UP};
    $tally{$tis}{DOWN}=0 if !$tally{$tis}{DOWN};
    my $t=$tally{$tis}{UP}+$tally{$tis}{DOWN};
    my $fr_AS=sprintf("%.2f",100*$t/$tallyTIS_AS{$tis});
    my $fr_UP=sprintf("%.2f",100*$tally{$tis}{UP}/$tallyTIS_AS{$tis});
    my $fr_DOWN=sprintf("%.2f",100*$tally{$tis}{DOWN}/$tallyTIS_AS{$tis});

    # no longer "-" for downregulated TS
    $data_print{$tis}="$tis\t$t\t$tally{$tis}{UP}\t$tally{$tis}{DOWN}\t$tallyTIS_AS{$tis}\t$fr_AS\t$fr_UP\t-$fr_DOWN\n";
    push(@order,"$fr_AS=$tis");
}
@order = sort {($b=~/(.+?)\=/)[0]<=>($a=~/(.+?)\=/)[0]}(@order);
#@order=sort{$b<=>$a}(@order);

foreach my $temp (@order){
    my ($tis)=$temp=~/\=(.+)/;
    print $data_print{$tis};
}


###### SUBROUTINES
sub parse_sample {
    my @t2=@_;
    my $ev=$t2[0];
    my $tis=$t2[1]; # actual tissuee
    my $AT=$t2[2]; # group
    my @t=@t2[3..$#t2];

    # quite convoluted, but maintained
    $sum{$ev}{$AT}+=$t[$val{$tis}];
    $tot{$ev}{$AT}++;
    if (defined $max{$ev}{$AT}){
	$max{$ev}{$AT}=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]>$max{$ev}{$AT};
    } else {
	$max{$ev}{$AT}=sprintf("%.2f",$t[$val{$tis}]);
    }
    if (defined $min{$ev}{$AT}){
	$min{$ev}{$AT}=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]<$min{$ev}{$AT};
    } else {
	$min{$ev}{$AT}=sprintf("%.2f",$t[$val{$tis}]);
    }
    $OK{$ev}{$AT}.="$tis,";
    $vals{$ev}{$AT}.="$t[$val{$tis}],";
}

sub median {
    my @t2=@_;
    my $ev=$t2[0];
    my $tis=$t2[1];

    my %V; my %med;
    @{$V{$tis}}=sort{$a<=>$b}(split(/\,/,$vals{$ev}{$tis}));
    my $ind="";
    $ind=$#{$V{$tis}}+1; 
    if ($ind%2!=0){
	my $m=int($ind/2);
	$med{$ev}{$tis}=$V{$tis}[$m];
    }
    else {
	my $m1=int($ind/2)-1;
	my $m2=int($ind/2);
	$med{$ev}{$tis}=($V{$tis}[$m1]+$V{$tis}[$m2])/2;
    }
    return $med{$ev}{$tis};
}

sub std_dev {
    my @t2=@_;
    my $ev=$t2[0];
    my $tis=$t2[1];
    my $av=$t2[2];
    my %SD;

    if (!$av){
	my $av=sprintf("%.2f",$sum{$ev}{$tis}/$tot{$ev}{$tis});
    }
    
    my @t3=split(/\,/,$vals{$ev}{$tis});
    my $sumCLAS=0;
    foreach my $n (@t3){
	$sumCLAS+=($n-$av)*($n-$av);
    }
    $SD{$ev}{$tis}=sprintf("%.2f",sqrt($sumCLAS/($#t3)));

    return $SD{$ev}{$tis};
}

sub std_dev2 {
    my @temp_vals=@_;

    my $sum1=0;
    foreach my $val (@temp_vals){
	$sum1+=$val;
    }
    my $av=$sum1/($#temp_vals+1);
    
    my $sumCLAS=0;
    foreach my $n (@temp_vals){
	$sumCLAS+=($n-$av)*($n-$av);
    }
    my $SD_others;
    $SD_others=sprintf("%.2f",sqrt($sumCLAS/($#temp_vals))) unless $#temp_vals==0;
    $SD_others="NA" if $#temp_vals==0;
    
    return $SD_others;
}

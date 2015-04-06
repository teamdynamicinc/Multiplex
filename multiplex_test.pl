#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::Primer3;
use Bio::Tools::Run::StandAloneBlastPlus;
                
#AUTHOR: Team Dynamic
#DATE STARTED: Mar 3/15
#DATE UPDATED: April 5/15
#PURPOSE: Produce primers from an inputted DNA sequence (fasta file) using primer3 tool
#NOTES: Must install BioPerl first following the instructions at http://www.bioperl.org/wiki/Installing_BioPerl_on_Windows
        #Must download Bio/Tools/Run and Bio/Roots/Roots modules and place in Perl>lib for primer3 tool to be accessed
        #Must copy primer3_core.exe in same location as THIS perl file to avoid "primer3 cannot be found" error
        #Use primer3 version 1.1.4 to avoid "missing SEQUENCE tag" error
        #Alter -path => in Primer3->new() to avoid "SH: Command not found" error
        #Copy primer3_config directory in same location as THIS perl file to avoid "thermodynamic approach" error
        #To avoid "can't locate object method 'new'" error, install correct version of Primer3.pm
        #Check output file (temp.out) to see actual left/right primer sequences

#*******************************************************************************************************#
#ATTACH JAVA GUI TO THIS PROGRAM FILE



#*******************************************************************************************************#
#USER-SELECTED ORGANISM FROM GUI USED AS INPUT



#*******************************************************************************************************#
#ACCESS ORGANISM'S FASTA FILE FROM DATABASE



#*******************************************************************************************************#
#PRIMER3 TOOL TO PRODUCE PRIMERS FROM FASTA FILE (Rebecca Allan, sourced from Chad Matsalla)
#Modified for multiplexing (Tiffany)
#Input a fasta file and create a new primer3 list

#results for primer3 will be put into "results" folder
#we do not want folder to have contents, to prevent duplications/confusions
opendir my $directory_result, "results" or die "can't open the directory\n"; #open directory
my @files_result = readdir $directory_result; #reads all file names in directory
closedir $directory_result;
        
foreach (@files_result){
    #if a file exists in directory
    if (-f "results/$_") { 
        die "your \"results\" folder has to be empty\n"; #will stop program
    }
}

my $x = 0; #counter for each result from primer3
my @results; #results from primer3

#java gui with fileupload will put fasta files from user into data folder
opendir my $directory_data_handler, "data" or die "can't open the directory\n"; #open directory
my @files_data = readdir $directory_data_handler; #reads all file names in directory
closedir $directory_data_handler;
            
my @fasta_data = grep{/\.fa$/} @files_data; #filters file names that ends only with .fa for fasta files
my @fasta_data_ordered = sort @fasta_data; #put in order
foreach (@fasta_data_ordered) {
                
    if (-f "data/$_") { #test if file exist is in Data directory
                                 
    my $seqio = Bio::SeqIO->new(-file=>"data/$_");
    my $seq = $seqio->next_seq;
    my @filename = split(/[\.\/]/,$_); #takes file extension and directory path off filename
    my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq,
                                                -outfile =>"results/temp_test_$filename[1].out", #output file named after fasta file name
                                                -path => "C:/Users/Tiffany/Desktop/Project/primer3_core"); #USE YOUR OWN PATH!
    
    #Test to see if primer3_core.exe is within the directory
    unless ($primer3->executable) {
        print STDERR "Primer3 can not be found. Is it installed?\n";
        exit(-1)
    }

#Adjust default values of specific arguments in primer3
$primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>65);
$primer3->add_targets('PRIMER_MIN_SIZE'=>20, 'PRIMER_MAX_SIZE'=>27, 'PRIMER_OPT_SIZE'=>20, 'PRIMER_DEFAULT_SIZE'=>20);
#Next line on primer product size range does not work, must be fixed
#$primer3->add_targets('PRIMER_PRODUCT_SIZE_RANGE'=>100..500);
#Number of primers produced
$primer3->add_targets('PRIMER_NUM_RETURN'=>5);

print "Please wait while ideal primer candidates being found in your sequence. This may take up to 5 minutes.\n";
$results[$x] = $primer3->run;

#Print the number of results (proves successful run)
#print "There were ", $results->number_of_results, " primers\n";

#*******************************************************************************************************#
#RETRIEVE PRIMERS FROM TEMP.OUT FILE (Rebecca Allan and Phuong Ma)
#Array to contain temp.out contents
my @file_contents;
#Arrays to contain the primers and their positions
my @left_primers = ();
my @left_primers_pos = ();
my @right_primers = ();
my @right_primers_pos = ();
#Array to contain amplimers for each primer pair
my @amplimers = ();

#Read the entire temp.out file into an array and then close file
open (my $filehandle, "<", "results/temp_test_$filename[1].out") or die "Error reading file.\n";

foreach (<$filehandle>) {
    push @file_contents, $_;
}
close $filehandle;
chomp @file_contents;

#Scalar value $temp_out for array @file_contents for regex
my $temp_out = join '', @file_contents;

for (my $i = 0; $i < $results[$x]->number_of_results; $i++) {
    #Perform different regex statements depending on the index value of each primer pair element
    if ($i > 0) {
        $temp_out =~ /(PRIMER_LEFT_$i\_SEQUENCE=)(.+)(PRIMER_RIGHT_$i\_SEQUENCE=)(.+)(PRIMER_LEFT_$i=)(.+)(\,.+)(PRIMER_RIGHT_$i=)(.+?)(\,.+)/;
        #Find primers and their positions, then push them into the correct array
        my $temp_left_seq = $2;
        my $temp_left_pos = $6;
        my $temp_right_seq = $4;
        my $temp_right_pos = $9;
        push @left_primers, $temp_left_seq;
        push @left_primers_pos, $temp_left_pos;
        push @right_primers, $temp_right_seq;
        push @right_primers_pos, $temp_right_pos;
        
        #Produce the amplimers and place in an array
        my $right = scalar reverse $temp_right_seq;
        $right =~ tr/AGCT/TCGA/;
        $temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    } else {
        $temp_out =~ /(PRIMER_LEFT_SEQUENCE=)(.+)(PRIMER_RIGHT_SEQUENCE=)(.+)(PRIMER_LEFT=)(.+)(\,.+)(PRIMER_RIGHT=)(.+?)(\,.+)/;
        #Find primers and their positions, then push them into the correct arra
        my $temp_left_seq = $2;
        my $temp_left_pos = $6;
        my $temp_right_seq = $4;
        my $temp_right_pos = $9;
        push @left_primers, $temp_left_seq;
        push @left_primers_pos, $temp_left_pos;
        push @right_primers, $temp_right_seq;
        push @right_primers_pos, $temp_right_pos;
        
        #Produce the amplimers and place in an array
        my $right = scalar reverse $temp_right_seq;
        $right =~ tr/AGCT/TCGA/;
        $temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    }
}

#Delete primers that have the same left/left or right/right primer positions
for(my $a = 0; $a < @left_primers; $a++) {
    for(my $b = 0; $b < @left_primers; $b++) {
        if ($a != $b) {
            if ($left_primers_pos[$a] == $left_primers_pos[$b]) {
                splice @left_primers, $b, 1;
                splice @left_primers_pos, $b, 1;
                splice @right_primers, $b, 1;
                splice @right_primers_pos, $b, 1;
                splice @amplimers, $b, 1;
            } elsif ($right_primers_pos[$a] == $right_primers_pos[$b]) {
                splice @left_primers, $b, 1;
                splice @left_primers_pos, $b, 1;
                splice @right_primers, $b, 1;
                splice @right_primers_pos, $b, 1;
                splice @amplimers, $b, 1;
            }
        }
    }
}

#open filehandler to append into txt file
open(my $fh_out2, ">>", "results/test.txt") or die "Could not open file $!";

#open output file to read
open (my $results_handler, "<", "results/temp_test_$filename[1].out") or die "Error reading file.\n";
while (my $results_output = <$results_handler>) {
    print $fh_out2 "$results_output"; #print all contents of file into appending filehandler for txt file
}

#print into textfile           
print $fh_out2 "Set: $filename[1] \n";
print "Set: $filename[1] \n";
#Display the primers, their position and the amplimers produced by them for the user
for (my $i = 0; $i < @left_primers; $i++) {                      
    printf $fh_out2 "%10s\t%25s\t%30s", "Primer Set $i:", "$left_primers[$i]", "$right_primers[$i]";
    printf "%10s\t%25s\t%30s", "Primer Set $i:", "$left_primers[$i]", "$right_primers[$i]";
    print $fh_out2 "\n";
    print "\n";
    printf $fh_out2 "%10s\t%25s\t%30s", "Position $i:", "$left_primers_pos[$i]", "$right_primers_pos[$i]";
    printf "%10s\t%25s\t%30s", "Position $i:", "$left_primers_pos[$i]", "$right_primers_pos[$i]";
    print $fh_out2 "\n";
    print "\n";
    printf $fh_out2 "%10s\t%25s", "Amplimer $i:", "$amplimers[$i]";
    printf "%10s\t%25s", "Amplimer $i:", "$amplimers[$i]";
    print $fh_out2 "\n \n";
    print "\n \n";    
}
close $fh_out2;
$x++; #counter increases for labelling results for next fasta file
}
}

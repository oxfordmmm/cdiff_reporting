#!/usr/bin/perl

#use DateTime;
use Time::Piece;

$raw = $ARGV[0];
$sample = "sample_name.tsv";

open(in, $raw);
open(sample, '>', $sample);

while(<in>){
    chomp;
    $f = $_;
    #print "$f\n";
    $date = localtime->strftime('%m/%d/%Y'); 
    #print "$datetime\n";
    @fname = split/\s/,$f;
    for $n (0..$#fname){}
    $fname[0] =~ s/\_R\_contigs\.fa//g;
    #print "Sample details\n";
    #print sample "Specimen identifier:\t$fname[0]\n";
    #print "Sequence GUUID:\t\t$fname[1]\t\tSequencing date:\t\t$date\n";
    #print sample "Report date:\t$date\n";
    #print sample "Multi-locus sequence type:\t$fname[2]\n";
    print "\nReport date\tSpecimen identifier\tMulti-locus sequence type\n\n";
    print "$date\t$fname[0]\t$fname[2]\n";
    print sample "Report date\tSpecimen identifier\tMulti-locus sequence type\n";
    print sample "$date\t$fname[0]\t$fname[2]\n";
}

close(in);
close(sample);

exit;

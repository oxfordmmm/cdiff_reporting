#!/usr/bin/perl

#use DateTime;
use Time::Piece;

$raw = $ARGV[0];
$sample = $ARGV[1];
$ribo = $ARGV[2];

open(in, $raw);
open(sample, '>', $sample);

while(<in>){
    chomp;
    $f = $_;
    #print "$f\n";
    $date = localtime->strftime('%d/%m/%Y'); 
    #print "$datetime\n";
    @fname = split/\s/,$f;
    for $n (0..$#fname){}
    $fname[0] =~ s/\_R\_contigs\.fa//g;
    #print "Sample details\n";
    #print sample "Specimen identifier:\t$fname[0]\n";
    #print "Sequence GUUID:\t\t$fname[1]\t\tSequencing date:\t\t$date\n";
    #print sample "Report date:\t$date\n";
    #print sample "Multi-locus sequence type:\t$fname[2]\n";
    #print "\nReport date\tSpecimen identifier\tMulti-locus sequence type\n\n";
    #print "$date\t$fname[0]\t$fname[2]\n";
}
open(ribo, $ribo);

print "Report date\tSpecimen identifier\tMLST\tRibotype\tCollection Date\tSource Hospital\n";
print sample "Report date\tSpecimen identifier\tMLST\tRibotype\tCollection Date\tSource Hospital\n";

while(<ribo>){
    chomp;
    $r = $_;
    #print "$r\n";
    @rt = split/\s/,$r;
    for $rt (0..$#rt){}
    #print "$rt[10]\n";
    #print "$fname[2]++\n";
    if ($fname[2] == $rt[0]){
        #print "$rt[0]\t$fname[2]\t$rt[1]\n";
        #print "$date\t$fname[0]\t$fname[2]\t$rt[10]\t\t\n";
        print "$date\t$fname[0]\t$fname[2]\t$rt[1]\t\t\n";
        print sample "$date\t$fname[0]\t$fname[2]\t$rt[1]\t\t\n";    
    }   
    else {
        #print "no match\n";
    }
}


close(ribo);
close(in);
close(sample);

exit;

#!/usr/bin/perl

#use PDF::API2;
#use DateTime;
use Time::Piece;

$fqc = $ARGV[0]; #FastQC output
$qst = $ARGV[1]; #Quast output

open(qc, $fqc);
open(qt, $qst);

$outqc = $ARGV[2];
#$sample = "sample_name.tsv";

#system("echo Sample Name\t$fqc | cut -f 1 -d '.' ");
#system("echo Sample Name\t$fqc | cut -f 1 -d '.' >> $outqc");

open(outqc, '>>', $outqc);
#open(sample, '>', $sample);

#print outqc "Sequence Quality\n";

print "\nGENOME QUALITY STATS\n";
print "\nMetric\tExpected value\tActual value\tQuality\n\n";
print outqc "Metric\tExpected value\tActual value\tQuality\n";

while(<qc>){
    chomp;
    $f = $_;
    #print "$f\n";
    #if ($f =~ "Basic Statistics" && $f =~ "pass"){print "Overall sequence quality\tPass\n"} elsif ($f =~ "Basic Statistics" && $f =~ "fail") {print "Overall sequencing quality\tFail\n"}
    if ($f =~ /Filename/){
        $date = localtime->strftime('%m/%d/%Y'); 
        #print "$datetime\n";
        @fname = split/\t/,$f;
        for $n (0..$#fname){}
        $fname[1] =~ s/\_raw\_reads\.fastq\.gz//g;
        #print sample "Sample details\n";
        #print sample "Specimen identifier:\t$fname[1]\tCollection date:\t$date\n";
        #print sample "Sequence GUUID:\t\t$fname[1]\tSequencing date:\t$date\n";
        #print sample "Sequence plate name:\t$fname[1]\ttReport date:\t$date\n";
    }
    if ($f =~ /Total Sequences/){
        @seqs = split/\t/,$f;
        for $s (0..$#seqs){}
        if ($seqs[1] == 1000000 | $seqs[1] < 20000000 | $seqs[1] == 20000000){
            $seqs[1] = sprintf("%.1f", $seqs[1]/1000000); 
            print "Total Sequences (M)\t1-20\t$seqs[1]\tPass\n";
            print outqc "Total Sequences (M)\t1-20\t$seqs[1]\tPass\n";
        }
        else {
            $seqs[1] = sprintf("%.1f", $seqs[1]/1000000); 
            print "Total Sequences (M)\t1-20\t$seqs[1]\tFail\n";
            print outqc "Total Sequences (M)\t1-20\t$seqs[1]\tFail\n";
        }
    }
    if ($f =~ /\%GC/){
        @gc = split/\t/,$f;
        for $g (0..$#gc){}
        if ($gc[1] == 27.9 | $gc[1] < 29.2 | $gc[1] == 29.2){
            print "%GC\t27.9-29.2\t$gc[1]\tPass\n";
            print outqc "%GC\t27.9-29.2\t$gc[1]\tPass\n";
        }
        else {
            print "%GC\t27.9-29.2\t$gc[1]\tFail\n";
            print outqc "%GC\t27.9-29.2\t$gc[1]\tFail\n";
        }
    }

    if ($f =~ /Sequence length/){
        @len = split/\t/,$f;
        # for $s (0..$#len){}
        #print "$len[1]*\n";
        # @seqlen = split/-/,$len[1];
        # for $s (0..$#seqlen){}
        # $mean = ($seqlen[0] + $seqlen[1])/2;
        #print "$mean**\n";
        $length = $len[1];
        if ($len[1] =~ /-/){
            @seqlen = split/-/,$len[1];
            $mean = ($seqlen[0] + $seqlen[1])/2;
            $length = $mean;
        }

        if ($length > 50 & $length < 150){
            print "Sequence length (bp)\t50-150\t$length\tPass\n";
            print outqc "Sequence length (bp)\t50-150\t$length\tPass\n";
        }
        else{
            print "Sequence length (bp)\t50-150\t$length\tFail\n";
            print outqc "Sequence length (bp)\t50-150\t$length\tFail\n";
        }
    }

}

close(qc);

while(<qt>){
    chomp;
    $q = $_;
    #print "$q\n";
    if ($q =~ "Total length" && $q =~ "(>= 0 bp)"){
        @gnom = split/\t/,$q;
        for $m (0..$#gnom){}
        if ($gnom[1] == 3900000 | $gnom[1] < 4500000 | $gnom[1] == 4600000){
            $gnom[1] = sprintf("%.1f", $gnom[1]/1000000); 
            print "Total assembly size (Mbp)\t3.9-4.5\t$gnom[1]\tPass\n";
            print outqc "Total assembly size (Mbp)\t3.9-4.5\t$gnom[1]\tPass\n";
        }
        else {
            $gnom[1] = sprintf("%.1f", $gnom[1]/1000000); 
            print "Total assembly size (Mbp)\t3.9-4.5\t$gnom[1]\tFail\n";
            print outqc "Total assembly size (Mbp)\t3.9-4.5\t$gnom[1]\tFail\n";
        }
    }
    if ($q =~ "Largest contig"){
        @con = split/\t/,$q;
        for $c (0..$#con){}
        if ($con[1] == 1000 | $con[1] < 1000000 | $con[1] == 1000000){
            $con[1] = sprintf("%.1f", $con[1]/1000);
            print  "Largest contig (Kbp)\t1-1000\t$con[1]\tPass\n";
            print outqc  "Largest contig (Kbp)\t1-1000\t$con[1]\tPass\n";
        }
        else {
            $con[1] = sprintf("%.1f", $con[1]/1000);
            print  "Largest contig (Kbp)\t1-1000\t$con[1]\tFail\n";
            print outqc  "Largest contig (Kbp)\t1-1000\t$con[1]\tFail\n";
        }   
    }

    if ($q =~ "N50"){
        @nfifty = split/\t/,$q;
        for $s (0..$#nfifty){}
        if ($nfifty[1] == 100000 | $nfifty[1] < 1000000 | $nfifty[1] == 1000000){
            $nfifty[1] = sprintf("%.1f", $nfifty[1]/1000);
            print "N50 (Kbp)\t10-1000\t$nfifty[1]\tPass\n";
            print outqc "N50 (Kbp)\t10-1000\t$nfifty[1]\tPass\n";
        }
        else {
            $nfifty[1] = sprintf("%.1f", $nfifty[1]/1000);
            print "N50 (Kbp)\t10-1000\t$nfifty[1]\tFail\n";
            print outqc "N50 (Kbp)\t10-1000\t$nfifty[1]\tFail\n";
        }
    }
    if ($q =~ /\%GC/){
        @gc = split/\t/,$q;
        for $g (0..$#gc){}
        if ($gc[1] == 27.9 && $gc[1] < 29.2 && $gc[1] == 29.2){
            print outqc "%GC\t27.9-29.2\t$gc[1]\tPass\n";
        }
        else {
            print outqc "%GC\t27.9-29.2\t$gc[1]\tFail\n";
        }
    }


}
print "\n\n";

#close(sample);

close(qt);

close(outqc);


exit;

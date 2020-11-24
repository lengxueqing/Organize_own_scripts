#!/usr/bin/perl -w
use warnings;
use strict;
use Data::Dumper;
##perl index_annovar.pl input.txt 1000
die "$0 <Annovar Database File> <BIN Size>" unless @ARGV ==2;
my $input_file =$ARGV[0];
my $bin_size =$ARGV[1];

if(! -e $input_file){
    die "$input_file not found!\n";
}


my $file_size= -s $input_file;

my %index;
open(my $in,"<",$input_file) or die "Could not open $input_file for indexing!\n";

my $previous_file_position = tell $in;
print $previous_file_position ."\n";

while(my $ln =<$in>){
    my($chr,$start,$stop)=split "\t",$ln;
    $chr=~s/^chr//i;
    my $bin_start = int($start/$bin_size)*$bin_size;
    #print $bin_start."\t";
    my $current_file_position =tell $in;
    #print "$current_file_position"."\n";
    
    if(!exists $index{$chr}->{$bin_start}){
        $index{$chr}->{$bin_start}=[$previous_file_position,$current_file_position];
    }else{
        $index{$chr}->{$bin_start}->[1]=$current_file_position;
    }
    $previous_file_position=$current_file_position;
 }

close $in;
print Dumper(\%index);

my $index_file=$input_file.".idx";
open(OUT,">$index_file");
print OUT "#BIN\t$bin_size\t$file_size\n";
foreach my $chr((1,10..19,2,20..22,3..9,"MT","X","Y")){
    foreach my $index_region(sort keys %{$index{$chr}}){
        my $start = $index{$chr}->{$index_region}->[0];
        my $stop = $index{$chr}->{$index_region}->[1];
        print OUT "$chr\t$index_region\t$start\t$stop\n";
    }
}
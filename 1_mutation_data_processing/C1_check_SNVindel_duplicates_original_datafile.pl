#!/usr/bin/perl -w 
#

open(FP,"../datasets/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf") or die; 
open(WR,">../datasets/PCAWG_filtered/final_consensus_passonly.snv_mnv_indel.icgc.public_duplicate_filt.maf") or die; 

$line=0; 
while($fp=<FP>) 
{
   chomp($fp); 
   $line++; 
   if($fp=~/^Hugo/) 
   {
     print WR "$fp\n";  
     next; 
   } 
   @arr=split(/\t/,$fp); 
   chop($arr[42]); 
   #$no=@arr;
   $id = $arr[0]."*$arr[1]*$arr[2]*$arr[3]*$arr[4]*$arr[7]*$arr[9]*$arr[42]";

   #print "$arr[42]\n"; 
   if(exists $list{$id}) 
   {
      print "$id $line $list{$id}\n"; 
   }
   else
   {
      print WR "$fp\n"; 
      $list{$id}=$line;
   }
   undef $id;
   undef(@arr); 
}
close(WR); 
close(FP); 

undef %list; 

#!/usr/bin/perl -w 
#
# Updated version 

## Read promoter information 
open(FP,"human_TFdata_programs/results/Homo_sapiens_GRCh37_Promoter_regions_1000.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Chromosome/) { next; } 
   @arr=split(/\s+/,$fp); 
   $chr=$arr[0]; 
   $id = $chr; #.$arr[5];
   if(exists $promlist{$id}) 
   {
      $cnt=$promlist{$id}[0];
      $promlist{$id}[1][$cnt]=$arr[6];
      $promlist{$id}[2][$cnt]=$arr[7];
      $promlist{$id}[3][$cnt]=$arr[1];
      $promlist{$id}[4][$cnt]=$arr[5]; ## strand
      $promlist{$id}[0]++; 
   }
   else
   {
      $promlist{$id}[1][0]=$arr[6]; 
      $promlist{$id}[2][0]=$arr[7];
      $promlist{$id}[3][0]=$arr[1]; 
      $promlist{$id}[4][0]=$arr[5]; ## strand
      $promlist{$id}[0]=1; 
   }
   undef $id; 
   undef(@arr); 
}
close(FP);

$intfile[0] = 'intronic_variants_Jung_etal2021/S2_intronic_vars_SS.txt';
$intfile[1] = 'intronic_variants_Jung_etal2021/S3_intronic_vars_proximal.txt';
$intfile[2] = 'intronic_variants_Jung_etal2021/S4_intronic_vars_deep.txt';
$intfile[3] = 'intronic_variants_Jung_etal2021/S5_intronic_vars_branchpoints.txt';
$intfile[4] = 'intronic_variants_Jung_etal2021/S6_exonic_vars.txt';

for($k=0;$k<5;$k++)
{
   $ncol = 5; 
   if($intfile[$k]=~/S2/) { $qno = 10; } 
   if($intfile[$k]=~/S3/) { $qno = 11; } 
   if($intfile[$k]=~/S4/) { $qno = 12; } 
   if($intfile[$k]=~/S5/) { $qno = 11; } 
   if($intfile[$k]=~/S6/) { $qno = 11; } 

   open(FP,$intfile[$k]) or die; 
   while($fp=<FP>) 
   {
      chomp($fp); 
      @arr=split(/\s+/,$fp); 
      $no=@arr;     

      if($fp=~/^chr/) 
      { 
	 for($li=$no-$ncol;$li<$no;$li++) 
	 {
	    #print "$arr[$li]:"; 
	    $header[$li]=$arr[$li]; 
	 }
	 undef(@arr);
	 undef $no; 
	 next; 
      } 

      $nval = "X"; 
      for($li=$no-5;$li<$no;$li++) 
      {
	 if($arr[$li] eq "O") 
	 {
	    $nval=$header[$li]; 
	    last; 
	 }
      } 

      #$qid = $arr[0]."*$arr[1]*$arr[$qno]*$arr[$qno-2]*$arr[$qno-1]"; 
      $qid = $arr[0]."*$arr[1]*$arr[$qno-2]*$arr[$qno-1]"; 
      $intlist{$qid}[0] = $nval; 
      $intlist{$qid}[1] = $arr[4]; 
      undef $qid; 

      undef $no; 
      undef(@arr); 
   }
   close(FP);
}

undef(@intfile);


open(FP,"../datasets/PCAWG_filtered/final_consensus_passonly.snv_mnv_indel.icgc.public_duplicate_filt.maf") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^Hugo\_Symbol/) 
   {
      next; 
   }
   chop($fp); 
   @arr=split(/\t/,$fp); 
   $chrlist{$arr[1]}=1; 
   $gene=$arr[0];
   $gene_prom = 'X'; 

   if(($arr[0] eq "Unknown" && $arr[5] eq "IGR") || ($arr[5] eq "5'Flank"))
   {
      $sid = $arr[1]; #.$arr[4];
      $sflag = 0;
      $gene = "Unknown"; 
      $gene_prom = 'X'; 
      $gene_strand = 'X'; 
      if($arr[5] eq "5'Flank") 
      {
	 $gene_prom = $arr[0]; 
      }
      if(exists $promlist{$sid}) 
      {
         for($lk=0;$lk<$promlist{$sid}[0];$lk++)
	 {
	    if($arr[2]>=$promlist{$sid}[1][$lk] && $arr[2]<=$promlist{$sid}[2][$lk])
	    {
	       $sflag = 1; 
	       $gene_prom = $promlist{$sid}[3][$lk]; 
	       $gene_strand = $promlist{$sid}[4][$lk]; 
	       last; 
	    }
	 }
      } 
      #if($sflag!=0)
      #{
	  #  print "$arr[1] $arr[2] $gene_prom $gene_strand $arr[42]\n"; #exit(); 
      #}
      undef $sid; 
   }

   ## Intronic variants comparison Jung et al., 2021 
   
   $intvar = 'X';
   $intgene = 'X'; 
   
   #$sid = $arr[1]."*$arr[2]*$arr[4]*$arr[7]*$arr[9]";  
   $sid = $arr[1]."*$arr[2]*$arr[7]*$arr[9]"; 

   if(exists $intlist{$sid}) 
   {
      $intvar = $intlist{$sid}[0]; 
      $intgene = $intlist{$sid}[1]; 
      print "$sid $intvar $intgene\n"; ##exit();  
   }
   undef $sid; 
   ## 


   $imp=$arr[5]; 
   $type=$arr[6]; 
   $ori=$arr[7];
   $mut=$arr[9];
   $can=$arr[41];
   $donid=$arr[42];
   #print "$can $donid $gene $imp $type $ori $mut\n";  

   if($type eq "SNP") 
   {
      ## $id=$can."*$donid*$ori*$mut"; 
      $id=$donid."*$ori*$mut"; 
      $muttyp{$ori."*$mut"}=1; 
      if(exists $mutlist{$id})
      {
         $cnt=$mutlist{$id}[0]; 
         $mutlist{$id}[1][$cnt]=$gene;    
         $mutlist{$id}[2][$cnt]=$gene_prom;    
         $mutlist{$id}[3][$cnt]=$imp;    
         $mutlist{$id}[4][$cnt]=$type;    
         $mutlist{$id}[5][$cnt]=$intvar;
         $mutlist{$id}[6][$cnt]=$intgene;
	 $mutlist{$id}[0]++; 
	 undef $cnt; 
      }
      else
      {
         $mutlist{$id}[1][0]=$gene;    
         $mutlist{$id}[2][0]=$gene_prom;    
         $mutlist{$id}[3][0]=$imp;    
         $mutlist{$id}[4][0]=$type;    
         $mutlist{$id}[5][0]=$intvar;
         $mutlist{$id}[6][0]=$intgene;
	 $mutlist{$id}[0]=1; 
      }
      undef $id; 

      if(exists $cntlist{$donid}) 
      {
	  $cntlist{$donid}++; 
      }
      else
      {
	  $cntlist{$donid}=1; 
      }
   }
   elsif($type eq "DNP") 
   {
      for($lk=0;$lk<length($ori);$lk++)
      {
	 $sub1=substr($ori,$lk,1); 
	 $sub2=substr($mut,$lk,1); 
      
	 ## $id=$can."*$donid*$sub1*$sub2"; 
	 $id=$donid."*$sub1*$sub2"; 

	 $muttyp{$sub1."*$sub2"}=1; 
         if(exists $mutlist{$id})
         {
            $cnt=$mutlist{$id}[0]; 
            $mutlist{$id}[1][$cnt]=$gene;    
            $mutlist{$id}[2][$cnt]=$gene_prom;    
            $mutlist{$id}[3][$cnt]=$imp;    
            $mutlist{$id}[4][$cnt]=$type;    
            $mutlist{$id}[5][$cnt]=$intvar;
            $mutlist{$id}[6][$cnt]=$intgene;
	    $mutlist{$id}[0]++; 
	    undef $cnt; 
         }
         else
         {
            $mutlist{$id}[1][0]=$gene;    
            $mutlist{$id}[2][0]=$gene_prom;    
            $mutlist{$id}[3][0]=$imp;    
            $mutlist{$id}[4][0]=$type;    
            $mutlist{$id}[5][0]=$intvar;
            $mutlist{$id}[6][0]=$intgene;
	    $mutlist{$id}[0]=1; 
         }
         undef $id; 

	 undef $sub1; 
	 undef $sub2; 

      }
      if(exists $cntlist{$donid}) 
      {
         $cntlist{$donid} += length($ori); 
      }
      else
      {
	  $cntlist{$donid} = length($ori); 
      }
   }
   elsif($type eq "INS")
   {
      $nn=length($mut); 
      for($lk=0;$lk<$nn;$lk++)
      {
	 $sub1="-"; 
	 $sub2=substr($mut,$lk,1); 
      
	 $muttyp{$sub1."*$sub2"}=1; 
	 
	 ## $id=$can."*$donid*$sub1*$sub2"; 
	 $id=$donid."*$sub1*$sub2"; 

         if(exists $mutlist{$id})
         {
            $cnt=$mutlist{$id}[0]; 
            $mutlist{$id}[1][$cnt]=$gene;    
            $mutlist{$id}[2][$cnt]=$gene_prom;    
            $mutlist{$id}[3][$cnt]=$imp;    
            $mutlist{$id}[4][$cnt]=$type;    
            $mutlist{$id}[5][$cnt]=$intvar;
            $mutlist{$id}[6][$cnt]=$intgene;
	    $mutlist{$id}[0]++; 
	    undef $cnt; 
         }
         else
         {
            $mutlist{$id}[1][0]=$gene;    
            $mutlist{$id}[2][0]=$gene_prom;    
            $mutlist{$id}[3][0]=$imp;    
            $mutlist{$id}[4][0]=$type;    
            $mutlist{$id}[5][0]=$intvar;
            $mutlist{$id}[6][0]=$intgene;
	    $mutlist{$id}[0]=1; 
         }
         undef $id; 

	 undef $sub1; 
	 undef $sub2; 
      }

      if(exists $cntlist{$donid}) 
      {
         $cntlist{$donid}+=$nn; 
      }
      else
      {
	 $cntlist{$donid}=$nn; 
      }
      undef $nn; 
   }
   elsif($type eq "DEL")
   {
      $nn=length($ori); 
      for($lk=0;$lk<$nn;$lk++)
      {
	 $sub1=substr($ori,$lk,1); 
	 $sub2="-"; 
      
	 $id=$donid."*$sub1*$sub2"; 
	 ## $id=$can."*$donid*$sub1*$sub2"; 
	 $muttyp{$sub1."*$sub2"}=1; 

         if(exists $mutlist{$id})
         {
            $cnt=$mutlist{$id}[0]; 
            $mutlist{$id}[1][$cnt]=$gene;    
            $mutlist{$id}[2][$cnt]=$gene_prom;    
            $mutlist{$id}[3][$cnt]=$imp;    
            $mutlist{$id}[4][$cnt]=$type;    
            $mutlist{$id}[5][$cnt]=$intvar;
            $mutlist{$id}[6][$cnt]=$intgene;
	    $mutlist{$id}[0]++; 
	    undef $cnt; 
         }
         else
         {
            $mutlist{$id}[1][0]=$gene;    
            $mutlist{$id}[2][0]=$gene_prom;    
            $mutlist{$id}[3][0]=$imp;    
            $mutlist{$id}[4][0]=$type;    
            $mutlist{$id}[5][0]=$intvar;
            $mutlist{$id}[6][0]=$intgene;
	    $mutlist{$id}[0]=1; 
         }
         undef $id; 

	 undef $sub1; 
	 undef $sub2; 
      }

      if(exists $cntlist{$donid}) 
      {
         $cntlist{$donid}+=$nn; 
      }
      else
      {
	  $cntlist{$donid}=$nn; 
      }
      undef $nn; 
   }

   undef $donid; 
   undef $can; 
   undef $mut; 
   undef $ori; 
   undef $type; 
   undef $imp; 
   undef $gene; 
   undef(@arr); 
}
close(FP); 


undef %intlist; 
undef %promlist; 


foreach $key (sort keys %chrlist) 
{
   print "$key\n"; 
}
undef %chrlist; 


foreach $key (sort keys %mutlist) 
{
   @qw=split(/\*/,$key);
   
   #$lid=$qw[0]."*$qw[1]"; 
   $lid = $qw[0]; 

   $pct=0; $cnn=0;  
   if(exists $cntlist{$lid}) 
   {
      $cnn=$cntlist{$lid};
      $pct=sprintf("%0.2f",$mutlist{$key}[0]/$cnn*100); 
   }
   #print "  $qw[1] $qw[2] $qw[3] $mutlist{$key}[0] $cnn $pct\n"; 

   if(exists $finlist{$lid}) 
   {
      $cnt=$finlist{$lid}[0];
      $finlist{$lid}[2][$cnt]=$qw[1];
      $finlist{$lid}[3][$cnt]=$qw[2];
      $finlist{$lid}[4][$cnt]=$pct;
      $finlist{$lid}[0]++;    
      undef $cnt; 
   }
   else
   {
      $finlist{$lid}[1]=$cnn; 
      $finlist{$lid}[2][0]=$qw[1];
      $finlist{$lid}[3][0]=$qw[2];
      $finlist{$lid}[4][0]=$pct;
      $finlist{$lid}[0]=1;    

      $finlist{$lid}[5]=0; ## Missense    
      $finlist{$lid}[6]=0; ## Nonsense    
      $finlist{$lid}[7]=0; ## Silent   
      $finlist{$lid}[8]=0; ## Intron 
      $finlist{$lid}[9]=0; ## Splice_site  
      $finlist{$lid}[10]=0; ## 3'UTR   
      $finlist{$lid}[11]=0; ## 5'UTR    
      #$finlist{$lid}[14]=0; ## 5'flank    
      $finlist{$lid}[12]=0; ## lincRNA     
      $finlist{$lid}[13]=0; ## inframe     
      $finlist{$lid}[14]=0; ## frameshift     
      $finlist{$lid}[15]=0; ## DeNovoInFrame
      $finlist{$lid}[16]=0; ## DeNovoOutofFrame
      $finlist{$lid}[17]=0; ## NumGenes    
      $finlist{$lid}[18]=0; ## NumGenes_StrEffMut    
      $finlist{$lid}[19]=0;  ## NumGenes_WeakEffMut 
      $finlist{$lid}[20]=0;  ## NumGenes_SilentMut
      $finlist{$lid}[21]=0;  ## NumGenes_TranscriptEffect
      $finlist{$lid}[22]=0;  ## NumProms
      $finlist{$lid}[23]=''; ## GenesAffected
      $finlist{$lid}[24]=''; ## Genes_StrEffMut
      $finlist{$lid}[25]=''; ## Genes_WeakEffMut
      $finlist{$lid}[26]=''; ## Genes_SilentMut
      $finlist{$lid}[27]=''; ## Genes_TranscriptEffMut
      $finlist{$lid}[28]=''; ## GenesPromAffected
      $finlist{$lid}[29]=0; ## NumFullIntRet
      $finlist{$lid}[30]=0; ## NumPartIntRet 
      $finlist{$lid}[31]=0; ## NumPseudoExonAct
      $finlist{$lid}[32]=0; ## NumFullExonSkip 
      $finlist{$lid}[33]=0; ## NumPartExonSkip 
      $finlist{$lid}[34]=''; ## GenesFullIntRet
      $finlist{$lid}[35]=''; ## GenesPartIntRet
      $finlist{$lid}[36]=''; ## GenesPseudoExonAct
      $finlist{$lid}[37]=''; ## GenesFullExonSkip
      $finlist{$lid}[38]=''; ## GenesPartExonSkip
   }

   $miss=0; $nons=0; $silent=0; $inframe=0; $frameshift=0; 
   $utr3pr=0; $utr5pr=0; $intron=0; $splsite=0;  
   $denoinframe=0; $denooutframe=0; $lincRNA = 0;  

   $numgene=0; $numgene_streff=0; $numgene_weakeff=0; $numgene_transeff=0; 
   $numprom=0; $numgene_silent=0;  

   $nfullintret = 0; $npartintret = 0; $npsdoexact=0; 
   $nfullexskip = 0; $npartexskip = 0; 

   
   for($lk=0;$lk<$mutlist{$key}[0];$lk++)
   {
      if((!exists $dup{$qw[0]."*$mutlist{$key}[1][$lk]"}) && ($mutlist{$key}[1][$lk] ne "Unknown"))
      {
         $dup{$qw[0]."*$mutlist{$key}[1][$lk]"}=1; 
	 #print "$mutlist{$key}[1][$lk]\t";
         $finlist{$lid}[23].="$mutlist{$key}[1][$lk],";    
	 $numgene++; 
      }
      if($mutlist{$key}[3][$lk] eq "Nonsense_Mutation" || $mutlist{$key}[3][$lk] eq "Frame_Shift_Del" || $mutlist{$key}[3][$lk] eq "Frame_Shift_Ins") 
      {
         if(!exists $sbdup1{$qw[0]."*$mutlist{$key}[1][$lk]"})
	 {
	    $sbdup1{$qw[0]."*$mutlist{$key}[1][$lk]"}=1; 
	    $numgene_streff++; 
	    $finlist{$lid}[24].="$mutlist{$key}[1][$lk],";
	 }
      } 
      if($mutlist{$key}[3][$lk] eq "Missense_Mutation" || $mutlist{$key}[3][$lk] eq "In_Frame_Del" || $mutlist{$key}[3][$lk] eq "In_Frame_Ins")
      {
         if(!exists $sbdup2{$qw[0]."*$mutlist{$key}[1][$lk]"})
	 {
	    $sbdup2{$qw[0]."*$mutlist{$key}[1][$lk]"}=1; 
	    $numgene_weakeff++; 
	    $finlist{$lid}[25].="$mutlist{$key}[1][$lk],";
	 }
      }
      if($mutlist{$key}[3][$lk] eq "Silent") 
      { 
         if(!exists $sbdup3{$qw[0]."*$mutlist{$key}[1][$lk]"})
	 {
	    $sbdup3{$qw[0]."*$mutlist{$key}[1][$lk]"}=1; 
	    $numgene_silent++; 
	    $finlist{$lid}[26].="$mutlist{$key}[1][$lk],"; 
	 }
      }
      if($mutlist{$key}[3][$lk] eq "Intron" || $mutlist{$key}[3][$lk] eq "Splice_Site" || $mutlist{$key}[3][$lk] eq "3'UTR" || $mutlist{$key}[3][$lk] eq "5'UTR")
      {
         if(!exists $sbdup4{$qw[0]."*$mutlist{$key}[1][$lk]"})
	 {
	    $sbdup4{$qw[0]."*$mutlist{$key}[1][$lk]"}=1; 
	    $numgene_transeff++;
	    $finlist{$lid}[27].="$mutlist{$key}[1][$lk],"; 
	 }
      }
      if((!exists $dup2{$qw[0]."*$mutlist{$key}[2][$lk]"}) && ($mutlist{$key}[2][$lk] ne "X"))
      {
         $dup2{$qw[0]."*$mutlist{$key}[2][$lk]"}=1; 
	 #print "$mutlist{$key}[2][$lk]\t";
         $finlist{$lid}[28].="$mutlist{$key}[2][$lk],";
         $numprom++; 	 
      }

      if($mutlist{$key}[5][$lk] eq "full_intron_ret") 
      {
	 $nfullintret++; 
         if(!exists $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"})
         {
	    $finlist{$lid}[34].="$mutlist{$key}[6][$lk],"; 
            $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"}=1
	 }
      }
      if($mutlist{$key}[5][$lk] eq "part_intron_ret") 
      {
         $npartintret++; 
         if(!exists $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"})
         {
            $finlist{$lid}[35].="$mutlist{$key}[6][$lk],"; 
            $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"}=1
	 }
      }
      if($mutlist{$key}[5][$lk] eq "pseudo_exon_act") 
      {
	 $npsdoexact++; 
         if(!exists $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"})
         {
	    $finlist{$lid}[36].="$mutlist{$key}[6][$lk],"; 
            $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"}=1
	 }
      }
      if($mutlist{$key}[5][$lk] eq "full_exon_skip") 
      {
	 $nfullexskip++;  
         if(!exists $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"})
         {
            $finlist{$lid}[37].="$mutlist{$key}[6][$lk],"; 
            $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"}=1
	 }
      }
      if($mutlist{$key}[5][$lk] eq "part_exon_skip") 
      {
	 $npartexskip++; 
         if(!exists $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"})
	 {
	    $finlist{$lid}[38].="$mutlist{$key}[6][$lk],"; 
            $dupint{$qw[0]."*$mutlist{$key}[5][$lk]*$mutlist{$key}[6][$lk]"}=1
	 }
      } 

      if($mutlist{$key}[3][$lk] eq "Missense_Mutation") { $miss++; } 
      if($mutlist{$key}[3][$lk] eq "Nonsense_Mutation") { $nons++; } 
      if($mutlist{$key}[3][$lk] eq "Silent") { $silent++; } 
      if($mutlist{$key}[3][$lk] eq "Intron") { $intron++; } 
      if($mutlist{$key}[3][$lk] eq "Splice_Site") { $splsite++; } 
      if($mutlist{$key}[3][$lk] eq "3'UTR") { $utr3pr++; } 
      if($mutlist{$key}[3][$lk] eq "5'UTR") { $utr5pr++; } 
      #if($mutlist{$key}[3][$lk] eq "5'Flank") { $flank5pr++; } 
      if($mutlist{$key}[3][$lk] eq "lincRNA") { $lincRNA++; } 
      if($mutlist{$key}[3][$lk] eq "Frame_Shift_Del" || $mutlist{$key}[3][$lk] eq "Frame_Shift_Ins") { $frameshift++; } 
      if($mutlist{$key}[3][$lk] eq "In_Frame_Del" || $mutlist{$key}[3][$lk] eq "In_Frame_Ins") { $inframe++; } 
      if($mutlist{$key}[3][$lk] eq "De_novo_Start_OutOfFrame") { $denooutframe++; } 
      if($mutlist{$key}[3][$lk] eq "De_novo_Start_InFrame") { $denoinframe++; } 
   }
   
   $finlist{$lid}[5]+=$miss; ## Missense    
   $finlist{$lid}[6]+=$nons; ## Nonsense    
   $finlist{$lid}[7]+=$silent; ## Silent   
   $finlist{$lid}[8]+=$intron; ## Intron 
   $finlist{$lid}[9]+=$splsite; ## Splice_site  
   $finlist{$lid}[10]+=$utr3pr; ## 3'UTR   
   $finlist{$lid}[11]+=$utr5pr; ## 5'UTR    
   $finlist{$lid}[12]+=$lincRNA; ## lincRNA     
   $finlist{$lid}[13]+=$inframe; ## inframe     
   $finlist{$lid}[14]+=$frameshift; ## frameshift     
   $finlist{$lid}[15]+=$denoinframe; ## DeNovoInFrame
   $finlist{$lid}[16]+=$denooutframe; ## DeNovoInFrame

   $finlist{$lid}[17]+=$numgene; ## NumGenes 
   $finlist{$lid}[18]+=$numgene_streff; ## NumGenes_StrEffMut
   $finlist{$lid}[19]+=$numgene_weakeff;  ## NumGenes_WeakEffMut
   $finlist{$lid}[20]+=$numgene_silent;  ## NumGenes_SilentMut 
   $finlist{$lid}[21]+=$numgene_transeff;  ## NumGenes_TranscriptEffect
   $finlist{$lid}[22]+=$numprom;  ## NumProms

   $finlist{$lid}[29]+=$nfullintret; ## NumFullIntRet
   $finlist{$lid}[30]+=$npartintret; ## NumPartIntRet 
   $finlist{$lid}[31]+=$npsdoexact; ## NumPseudoExonAct
   $finlist{$lid}[32]+=$nfullexskip; ## NumFullExonSkip 
   $finlist{$lid}[33]+=$npartexskip; ## NumPartExonSkip 

   undef $cnn; 
   undef $pct; 
   undef $lid; 
   undef(@qw);  
}

undef %dupint; 
undef %dup; 
undef %dup2; 
undef %sbdup1; 
undef %sbdup2; 
undef %sbdup3; 
undef %sbdup4; 
undef %cntlist; 
undef %mutlist; 

open(WR,">../results_June2025/MutationData_DonorWise.txt") or die; 
open(WR2,">../results_June2025/MutationGeneList_DonorWise.txt") or die; 
print WR "Cancer\tDonorID\t"; 
print WR2 "Cancer\tDonorID\t"; 
print WR "TumorType_WGS  project_code  donor_sex  donor_vital_status first_therapy_type first_therapy_response donor_age_at_diagnosis_years donor_survival_time_days  tumour_stage  MT_CopyNumber MT_SNV  MT_INDEL Tumor_Purity  Tumor_Ploidy NUC_SNV NUC_INDEL Ratio_MTCopyNumber_CancerToNormal\t"; 
print WR "Total_Num_Mutations\tTot_Nuc_Mut\t"; 

foreach $lty (sort keys %muttyp)
{
   @qw=split(/\*/,$lty);
   print WR "Perc_$qw[0]>$qw[1]\t"; 
   undef(@qw); 
}
print WR "\tPercSNP  PercINS  PercDEL\t"; 
#print WR "PercMissense PercNonsense PercSilent PercIntron PercSpliceSite  Perc3'UTR  Perc5'UTR PerclincRNA PercinFrameIndel PercframeShiftIndel PercdenovoStartInframe PercdenovoStartOutframe\t"; 
print WR "NumMissense NumNonsense NumSilent NumIntron NumSpliceSite  Num3'UTR  Num5'UTR NumlincRNA NuminFrameIndel NumframeShiftIndel NumdenovoStartInframe NumdenovoStartOutframe\t"; 
print WR "NumGenesAffected\tNumGenes_StrMutEff\tNumGenes_WeakMutEff\tNumGenes_SilentMut\tNumGenes_TranscrEff\tNumPromsAffected\tNumGenesPromsAffected\t";
print WR "NumFullIntRet\tNumPartIntRet\tNumPsdoExonAct\tNumFullExonSkip\tNumPartExonSkip\n";
print WR2 "GenesAffected\tGenesStrMutEff\tGenesWeakMutEff\tGenesSilentMut\tGenesTranscrEff\tPromotersAffected\tGenesCdsPromsAffected\t"; 
print WR2 "GenesFullIntRet\tGenesPartIntRet\tGenesPsdoExonAct\tGenesFullExonSkip\tGenesPartExonSkip\n";

open(FP,"../results_June2025/PCAWG_MitoData_collected.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\s+/,$fp); 
   if($fp=~/^TumorType/) 
   {
	   ##print "$arr[0] $arr[2] $arr[3] $arr[6] $arr[21] $arr[25] $arr[28] $arr[29] $arr[40]\n"; 
      undef(@arr); 
      next; 
   } 
   $mtdatalist{$arr[3]}[0]=$arr[0];
   $mtdatalist{$arr[3]}[1]=$arr[2];
   $mtdatalist{$arr[3]}[2]=$arr[6];
   $mtdatalist{$arr[3]}[3]=$arr[9];
   $mtdatalist{$arr[3]}[4]=$arr[10];
   $mtdatalist{$arr[3]}[5]=$arr[11];
   $mtdatalist{$arr[3]}[6]=$arr[21];
   $mtdatalist{$arr[3]}[7]=$arr[25];
   $mtdatalist{$arr[3]}[8]=$arr[26];
   $mtdatalist{$arr[3]}[9]=$arr[27];
   $mtdatalist{$arr[3]}[10]=$arr[28];
   $mtdatalist{$arr[3]}[11]=$arr[29];
   $mtdatalist{$arr[3]}[12]=$arr[31];
   $mtdatalist{$arr[3]}[13]=$arr[32];
   $mtdatalist{$arr[3]}[14]=$arr[40];

   $mtdatalist{$arr[3]}[15]=$arr[7]; ## donor_vital_status
   $mtdatalist{$arr[3]}[16]=$arr[12]; ## donor_survival_time_days
   $mtdatalist{$arr[3]}[17]=$arr[41]; ## Total_number_mutations
   $mtdatalist{$arr[3]}[18]=$arr[3]; ## DonorID
   undef(@arr); 
}
close(FP);

foreach $key (sort keys %finlist)
{
   #@qw=split(/\*/,$key); 
   #print WR "$qw[0]\t$qw[1]\t";
   #print WR2 "$qw[0]\t$qw[1]\t";

   if(exists $mtdatalist{$key})
   {
      print WR "$mtdatalist{$key}[1] $mtdatalist{$key}[18] $mtdatalist{$key}[0]  $mtdatalist{$key}[1]  $mtdatalist{$key}[2]  $mtdatalist{$key}[15]  $mtdatalist{$key}[3]  $mtdatalist{$key}[4]  $mtdatalist{$key}[5]  $mtdatalist{$key}[16]  $mtdatalist{$key}[6]  $mtdatalist{$key}[7]  $mtdatalist{$key}[8]  $mtdatalist{$key}[9]  $mtdatalist{$key}[10]  $mtdatalist{$key}[11] $mtdatalist{$key}[12] $mtdatalist{$key}[13] $mtdatalist{$key}[14] $mtdatalist{$key}[17]\t";
   }
   else { next; } 

   print WR "$finlist{$key}[1]\t";

   foreach $lty (keys %muttyp) 
   {
      $muttyp{$lty}=0; 
   }
   undef $lty; 

   for($lk=0;$lk<$finlist{$key}[0];$lk++)
   {
      $loc=$finlist{$key}[2][$lk]."*$finlist{$key}[3][$lk]"; 
      $muttyp{$loc}=$finlist{$key}[4][$lk]; 
      #print "$finlist{$key}[2][$lk]  $finlist{$key}[3][$lk]  $finlist{$key}[4][$lk]\t"; 
   }

   $ins=0; $snp=0; $del=0; 
   
   foreach $lty (sort keys %muttyp)
   {
      @we=split(/\*/,$lty);
      print WR "$muttyp{$lty}\t"; 
      if($we[0] eq "-") { $ins+=$muttyp{$lty}; }
      elsif($we[1] eq "-") { $del+=$muttyp{$lty}; }
      else { $snp+=$muttyp{$lty}; }
      undef(@we);
   }
   undef $lty; 

   print WR "\t$snp  $ins  $del\t";

   for($li=5;$li<=16;$li++) 
   {
      print WR "$finlist{$key}[$li]\t";
      #$perc = sprintf("%0.4f",$finlist{$key}[$li]/$finlist{$key}[1]*100); 
      #print WR "$perc\t"; 
   }
   for($li=17;$li<=22;$li++) 
   {
      print WR "$finlist{$key}[$li]\t"; 
   }

   $genepromlist = $finlist{$key}[23]; 
   $numtot = 0; 
   @tr = split(/\,/,$finlist{$key}[23]);
   $numgene = @tr; 
   for($lk=0;$lk<$numgene;$lk++) 
   {
       $duplist{$tr[$lk]}=1;
   }
   $numtot=$numgene; 	
   undef(@tr);

   @tr = split(/\,/,$finlist{$key}[28]);
   $numprom = @tr; 
   for($lk=0;$lk<$numprom;$lk++) 
   {
      if(!exists $duplist{$tr[$lk]})
      {       
         $duplist{$tr[$lk]}=1;
	 $genepromlist.="$tr[$lk],"; 
         $numtot++; 	
	 
      }
   }
   undef(@tr);
   print WR "$numtot\t";


   for($li=29;$li<=33;$li++) 
   {
      print WR "$finlist{$key}[$li]\t"; 
   }
   print WR "\n"; 


   if($finlist{$key}[23] eq '') { $finlist{$key}[23]='X'; } 
   if($finlist{$key}[24] eq '') { $finlist{$key}[24]='X'; } 
   if($finlist{$key}[25] eq '') { $finlist{$key}[25]='X'; } 
   if($finlist{$key}[26] eq '') { $finlist{$key}[26]='X'; } 
   if($finlist{$key}[27] eq '') { $finlist{$key}[27]='X'; } 
   if($finlist{$key}[28] eq '') { $finlist{$key}[28]='X'; } 

   print WR2 "$finlist{$key}[23]\t$finlist{$key}[24]\t$finlist{$key}[25]\t$finlist{$key}[26]\t$finlist{$key}[27]\t$finlist{$key}[28]\t$genepromlist\t";

   if($finlist{$key}[34] eq '') { $finlist{$key}[34]='X'; } 
   if($finlist{$key}[35] eq '') { $finlist{$key}[35]='X'; } 
   if($finlist{$key}[36] eq '') { $finlist{$key}[36]='X'; } 
   if($finlist{$key}[37] eq '') { $finlist{$key}[37]='X'; } 
   if($finlist{$key}[38] eq '') { $finlist{$key}[38]='X'; } 

   print WR2 "$finlist{$key}[34]\t$finlist{$key}[35]\t$finlist{$key}[36]\t$finlist{$key}[37]\t$finlist{$key}[38]\n";

   undef $numprom;
   undef $numgene;
   undef $genepromlist; 
   #undef(@qw); 
}
close(WR2);
close(WR); 

undef %duplist; 
undef %mtdatalist; 
undef %finlist; 
undef %muttyp; 

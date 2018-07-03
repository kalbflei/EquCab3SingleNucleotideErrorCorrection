/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.louisville.cgemm.twilightassembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 *
 * @author tedkalbfleisch
 */
public class FinalAlleleReplacement {
    
    public static void main(String[] args){
        
        try{
            HashMap TwilightSNPs = new HashMap();
            TwilightSNPs = getTwilightSNPs(TwilightSNPs,"/scratch/large/tskalb01/genomes/Equus_caballus/Ec_build-3k/definietlyFix.vcf");
            TwilightSNPs = getTwilightSNPs(TwilightSNPs,"/scratch/large/tskalb01/genomes/Equus_caballus/Ec_build-3k/probablyFix.vcf");
                        
            String line = null;
            StringTokenizer st = null;
            String chromosomeName = null;
            String position = null;
            
            HashMap outputHash = new HashMap();
            HashMap errHash    = new HashMap();
            HashMap contigHash = getContigHash();
            
            StringBuffer key = null;
            
            BufferedReader bufferedReader = new BufferedReader(new FileReader(new File("Ec_build-3_k.update.fasta")));
            boolean readFlag = false;
            StringBuffer seqBuffer = null;
            String seqName = null;
            while((line=bufferedReader.readLine())!=null){
                if(line.startsWith(">")){
                    chromosomeName = line.replace(">", "").trim();
                    if(contigHash.containsKey(seqName)){
                        processSeq(seqName,seqBuffer.toString(),TwilightSNPs);
                    }else{
                        if(seqBuffer!=null){
                            System.out.println(">" + seqName);
                            System.out.println(seqBuffer.toString());
                            System.out.flush();
                        }
                    }
                    seqBuffer = new StringBuffer();
                }else{
                    seqName = chromosomeName;
                    seqBuffer.append(line);
                }
                
            }
            if(contigHash.containsKey(seqName)){
                processSeq(seqName,seqBuffer.toString(),TwilightSNPs);
            }else{
                System.out.println(">" + seqName);
                System.out.println(seqBuffer.toString());
                System.out.flush();
            }
            
            bufferedReader.close();
            
        }catch(IOException ioException){
           ioException.printStackTrace();
        }
        
    }
    
    private static HashMap getContigHash(){
        
        //The only sequences that were modified were those that corresponded to the incorporated chromosomes
        //chr1-chrX
        
        HashMap contigHash = new HashMap();
        contigHash.put("Contig3799",null);
        contigHash.put("Contig3804",null);
        contigHash.put("Contig3805",null);
        contigHash.put("Contig4195",null);
        contigHash.put("Contig4198",null);
        contigHash.put("Contig4664",null);
        contigHash.put("Contig4301",null);
        contigHash.put("Contig3800",null);
        contigHash.put("Contig3803",null);
        contigHash.put("Contig4085",null);
        contigHash.put("Contig3802",null);
        contigHash.put("Contig4197",null);
        contigHash.put("Contig3801",null);
        contigHash.put("Contig4194",null);
        contigHash.put("Contig4202",null);
        contigHash.put("Contig4199",null);
        contigHash.put("Contig4196",null);
        contigHash.put("Contig4451",null);
        contigHash.put("Contig4300",null);
        contigHash.put("Contig4388",null);
        contigHash.put("Contig4302",null);
        contigHash.put("Contig4303",null);
        contigHash.put("Contig4201",null);
        contigHash.put("Contig4304",null);
        contigHash.put("Contig4200",null);
        contigHash.put("Contig4299",null);
        contigHash.put("Contig4450",null);
        contigHash.put("Contig4567",null);
        contigHash.put("Contig4387",null);
        contigHash.put("Contig4449",null);
        contigHash.put("Contig4568",null);
        contigHash.put("Contig4759",null);
        contigHash.put("Contig3798",null);
              
        return contigHash;
        
    }
    
    public static void processSeq(String seqName,String seqBuffer,HashMap TwilightSNPs){
        
        StringBuffer seqOut = new StringBuffer();
        int position = 0;
        String allele = null;
        String refAllele = null;
        int delta = 0;
        int val = 0;
        int counter = 0;
        int diffCounter = 0;
        for(int i=0;i<seqBuffer.length();){
            position = i+1;
            if((TwilightSNPs.containsKey(seqName + ":" + position))){
                

                allele = getVarAllele((String)TwilightSNPs.get(seqName + ":" + position));
                refAllele = getRefAllele((String)TwilightSNPs.get(seqName + ":" + position));
                i+=refAllele.length();
                seqOut.append(allele);
                if(allele.length()!=refAllele.length()){
                    delta+= allele.length()-refAllele.length();
                    System.err.println(delta);
                }
                diffCounter++;
            }else{
  
                seqOut.append(seqBuffer.charAt(i));
                i++;
            }

        }
        
        System.out.println(">" + seqName);
        System.out.println(seqOut.toString());
        System.out.flush();
        System.err.println("!!!!" + delta);
        System.err.println("diffCount=" + diffCounter);
    }
    
    public static String getVarAllele(String vcfRecord){
        
        StringTokenizer st = new StringTokenizer(vcfRecord,"\t,");
        st.nextToken();
        st.nextToken();
        st.nextToken();
        String[] alleles = new String[2];
        alleles[0] = st.nextToken();
        alleles[1] = st.nextToken();
        String token = null;
        /*
        while(st.hasMoreTokens()){
            token = st.nextToken();
        }
        System.err.println("Token=" + token);
        System.err.flush();
        int index = Integer.parseInt(token);
        
        if(index>1){
            System.err.println(vcfRecord);
            return alleles[0];
        }
        */
        return alleles[1];
         
    }
    
    public static String getRefAllele(String vcfRecord){
        
        StringTokenizer st = new StringTokenizer(vcfRecord,"\t,");
        String chromosomeName = st.nextToken();
        String position = st.nextToken();
        String varName = st.nextToken();
        String refAllele = st.nextToken();

        return refAllele;
         
    }
    
    public static String getErrAllele(String vcfRecord){
        
        StringTokenizer st = new StringTokenizer(vcfRecord,"\t");
        String token = null;
        st.nextToken();
        st.nextToken();
        st.nextToken();
        String refAllele = st.nextToken();
        String varAllele = st.nextToken();
        for(int i=0;i<5;i++){
            
            token = st.nextToken();
        }
        
        st = new StringTokenizer(token,":");
        token = st.nextToken();
        if(token.equals("1/1") || token.equals("1|1")){
            st = new StringTokenizer(varAllele,",");
            return st.nextToken();
        }else{
            return refAllele;
        }

    }
    
    private static HashMap getTwilightSNPs(HashMap twilightHash, String filename) throws IOException{
        
        BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(filename)));
                
        String line = null;
        StringTokenizer st = null;
        String chromosomeName = null;
        String position = null;
        String token = null;
        
        while((line=bufferedReader.readLine())!=null){
            if(line.startsWith("#")){
                continue;
            }
            
            st = new StringTokenizer(line);
            chromosomeName = st.nextToken();
            position = st.nextToken();
            
            for(int i=0;i<8;i++){
                token = st.nextToken();
            }
            
            //st = new StringTokenizer(token,":");
            if(!twilightHash.containsKey(chromosomeName + ":" + position)){
                twilightHash.put(chromosomeName + ":" + position, line);
            }else{
                System.err.println(line);
                System.err.flush();
                System.exit(1);
            }
            //System.out.println(token);
        }
        System.out.flush();

        return twilightHash;
        
    }
    
    private static boolean compareAlleles(String vcfRecord,String twilightGenotypeString){
        
        StringTokenizer st = new StringTokenizer(vcfRecord);
        String token = null;
        for(int i=0;i<10;i++){
            token = st.nextToken();
        }
        
        st = new StringTokenizer(token,":");
        String tenXallele = st.nextToken();
        st = new StringTokenizer(twilightGenotypeString,":/");
        
        for(int i=0;i<2;i++){
            if(!tenXallele.contains(st.nextToken())){
                System.err.println(vcfRecord);
                System.err.flush();
                return false;
            }
        }
        
        return true;
    }
    
}

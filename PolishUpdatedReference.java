/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.louisville.cgemm.twilightassembly;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 *
 * @author tedkalbfleisch
 */
public class PolishUpdatedReference {
    
    public static void main(String[] args){
        
        try{
            
            //HashMap TwilightSNPs = getTwilightSNPs("/scratch/large/tskalb01/Twilight/MacLeod_PCR_Free/WGS/MacLeod_PCR_Free.Ec_build-3k.vcf");
            
            BufferedReader bufferedReader = new BufferedReader(new FileReader(new File("thoroughbredVsTwilight.Ec_build-3k_update.vcf")));
            BufferedWriter definitelyFix = new BufferedWriter(new FileWriter(new File("definietlyFix.vcf")));
            BufferedWriter possiblyFix   = new BufferedWriter(new FileWriter(new File("probablyFix.vcf")));
            BufferedWriter remainingProblems   = new BufferedWriter(new FileWriter(new File("RemainingProblems.vcf")));
            String line = null;
            StringTokenizer st = null;
            String chromosomeName = null;
            String position = null;
            
            HashMap outputHash = new HashMap();
            HashMap errHash    = new HashMap();
            HashMap contigHash = getContigHash();
            String[] animalNames = null;
            String[] genotypeString = null;
            StringBuffer key = null;
            boolean homozygoteStrict = false;
            boolean homozygoteLenient = false;
            String genotype;
            
            while((line=bufferedReader.readLine())!=null){
                
                if(line.startsWith("#CHROM")){
                    st = new StringTokenizer(line);
                    animalNames = new String[(st.countTokens()-9)];
                    for(int i=0;i<9;i++){
                        st.nextToken();
                    }
                    for(int i=0;i<animalNames.length;i++){
                        animalNames[i]=st.nextToken();
                        System.out.println(animalNames[i]);
                    }                    
                }
                if(line.startsWith("#")){
                    definitelyFix.write(line);
                    definitelyFix.newLine();
                    possiblyFix.write(line);
                    possiblyFix.newLine(); 
                    remainingProblems.write(line);
                    remainingProblems.newLine();
                            
                    
                    continue;
                }

                st = new StringTokenizer(line);
                
                if(st.countTokens()<=3){
                    continue;
                }
                
                genotypeString = new String[(st.countTokens()-9)];
                chromosomeName = st.nextToken();
                if(!contigHash.containsKey(chromosomeName)){
                    continue;
                }
                for(int i=1;i<9;i++){
                    st.nextToken();
                }
                for(int i=0;i<genotypeString.length;i++){
                    genotypeString[i]=st.nextToken();                   
                }                 
                
                homozygoteStrict = true;
                homozygoteLenient = true;
                
                for(int i=0;i<genotypeString.length;i++){
                    st = new StringTokenizer(genotypeString[i],":");
                    genotype = st.nextToken();
                    if(!genotype.equals("1/1")){
                        homozygoteStrict = false;
                    }
                    if(!(genotype.equals("./.") || genotype.equals("1/1"))){
                        homozygoteLenient = false;
                    }
                }
                
                if(homozygoteStrict){
                    definitelyFix.write(line);
                    definitelyFix.newLine();
                }else if(homozygoteLenient && !homozygoteStrict){                   
                    possiblyFix.write(line);
                    possiblyFix.newLine(); 
                }else if(genotypeString[genotypeString.length -1].equals("1/1") && !homozygoteLenient && !homozygoteStrict){
                    remainingProblems.write(line);
                    remainingProblems.newLine();
                }
                
                
            }
            
            definitelyFix.close();
            possiblyFix.close();
            remainingProblems.close();
            
        }catch(IOException ioException){
           ioException.printStackTrace();
        }
        
    }
    
    private static HashMap getContigHash(){
        
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
    
    public static void processSeq(String seqName,String seqBuffer,HashMap TwilightSNPs, HashMap outputHash,HashMap errHash){
        
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
            if((outputHash.containsKey(seqName + ":" + position) && TwilightSNPs.containsKey(seqName + ":" + position))){
                
                if(compareAlleles((String)outputHash.get(seqName + ":" + position),(String)TwilightSNPs.get(seqName + ":" + position))){

                    if(outputHash.containsKey(seqName + ":" + position)){
                        allele = getAllele((String)outputHash.get(seqName + ":" + position));
                        refAllele = getRefAllele((String)outputHash.get(seqName + ":" + position));
                        i+=refAllele.length();
                        seqOut.append(allele);
                        if(allele.length()!=refAllele.length()){
                            delta+= allele.length()-refAllele.length();
                            System.err.println(delta);
                        }

                    }else if(errHash.containsKey(seqName + ":" + position)){
                        allele = getErrAllele((String)errHash.get(seqName + ":" + position));
                        refAllele = getRefAllele((String)errHash.get(seqName + ":" + position));
                        i+=refAllele.length();
                        seqOut.append(allele);
                        if(allele.length()!=refAllele.length()){
                            delta+= allele.length()-refAllele.length();
                            System.err.println(delta);
                        }

                    }else{
                        seqOut.append(seqBuffer.charAt(i));
                        i++;
                    }
                    
                }else{
                    diffCounter++;
                    seqOut.append(seqBuffer.charAt(i));
                    i++;
                }
                
            }else{
                if((outputHash.containsKey(seqName + ":" + position) || TwilightSNPs.containsKey(seqName + ":" + position)))
                    diffCounter++;
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
    
    public static String getAllele(String vcfRecord){
        
        StringTokenizer st = new StringTokenizer(vcfRecord,"\t,");
        st.nextToken();
        st.nextToken();
        st.nextToken();
        String[] alleles = new String[2];
        alleles[0] = st.nextToken();
        alleles[1] = st.nextToken();
        String token = null;
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
        
        return alleles[index];
         
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

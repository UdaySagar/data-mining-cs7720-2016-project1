/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package Panjala_Project1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;


/**
 *
 * @author Uday Sagar Panjala - U00771807
 */
public class CorrelationCoeff {
    HashMap<Integer, Double> hash = new HashMap<>();
     ArrayList lines = new ArrayList(); 
     static String[] fields;
     String[][] strings;
     Double correlation_coeff;

    public CorrelationCoeff(String data_Path,int k, ArrayList<Integer> gene_set1, ArrayList<Integer> gene_set2) throws FileNotFoundException, IOException
    {
      BufferedWriter out_cc = new BufferedWriter(new FileWriter(System.getProperty("user.dir") + "/src/Panjala_Project1/Output_CorrelationCoefficient.txt"));                     
        BufferedReader br=new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(data_Path)));
    	out_cc.write( "################################################################\n" );          
    	out_cc.write( "Correlation Coefficient of the selected Top K Genes" +"\n" );
    	out_cc.write( "################################################################\n" );  
    	System.out.println( "###############################################################" );                
    	System.out.println( "Correlation Coefficient of the selected Top K Genes" +"\n" );        
    	System.out.println( "###############################################################" );          
        for(String line = br.readLine();line != null;line = br.readLine()) {
            
                fields = line.split(", ");
            
                lines.add(fields);
             }
         strings= (String[][]) lines.toArray(new String[lines.size()][]);
                 int i,j,q;

         for( i=0;i<k;i++){
              for( q=0;q<k;q++){             
                 double[] gene1_values = new double[strings.length];
                 double[] gene2_values = new double[strings.length];

//////// To fill first gene all instance values                 
            for(j=0;j<strings.length;j++){
            if(strings[j][gene_set1.get(i)].contains(",")){
                    String data_class[]=strings[j][gene_set1.get(i)].split(",");
                    hash.put(j,Double.parseDouble(data_class[0]));
                    double val = Double.parseDouble(data_class[0]);  
                    gene1_values[j] = val;
                    
                }
            else{
                    double val = Double.parseDouble(strings[j][gene_set1.get(i)]);  
                    gene1_values[j] = val;               
            }
            }
//////// To fill second gene all instance values
            for(j=0;j<strings.length;j++){
            if(strings[j][gene_set2.get(q)].contains(",")){
                    String data_class[]=strings[j][gene_set2.get(q)].split(",");
                    hash.put(j,Double.parseDouble(data_class[0]));
                    double val = Double.parseDouble(data_class[0]);  
                    gene2_values[j] = val;
                    
                }
            else{
                    double val = Double.parseDouble(strings[j][gene_set2.get(q)]);  
                    gene2_values[j] = val;               
            }
            }
/////////// Done with filling /////////////////            
      Output2File("Correlation Coefficient of Gene no. "+gene_set1.get(i)+" and Gene no. "+gene_set2.get(q)+" is "+findCorrelationCoefficient(gene1_values, gene2_values)+"\n", out_cc);
              }        
}
         


        out_cc.close();
         
    }

  public static double findCorrelationCoefficient(double[] gene1_set_doublevalues, double[] gene2_set_doublevalues) {

    double standardDevOfX = 0.0;
    double standardDevOfY = 0.0;
    double standardDevOfXX = 0.0;
    double standardDevOfYY = 0.0;
    double standardDevOfXY = 0.0;

    int no_of_geneInstances = gene1_set_doublevalues.length;

    for(int i = 0; i < no_of_geneInstances; ++i) {
      double x = gene1_set_doublevalues[i];
      double y = gene2_set_doublevalues[i];

      standardDevOfX += x;
      standardDevOfY += y;
      standardDevOfXX += x * x;
      standardDevOfYY += y * y;
      standardDevOfXY += x * y;
    }

    // Finding covariation
    double covariation = standardDevOfXY / no_of_geneInstances - standardDevOfX * standardDevOfY / no_of_geneInstances / no_of_geneInstances;
    // standard error of x
    double X_Sigma = Math.sqrt(standardDevOfXX / no_of_geneInstances -  standardDevOfX * standardDevOfX / no_of_geneInstances / no_of_geneInstances);
    // Finding standard error
    double Y_Sigma = Math.sqrt(standardDevOfYY / no_of_geneInstances -  standardDevOfY * standardDevOfY / no_of_geneInstances / no_of_geneInstances);

    // Finding Correlation Coefficient
    return covariation / X_Sigma / Y_Sigma;
  }    
    
  public void Output2File(String write_Text, BufferedWriter out_cc){
         try{
           out_cc.write(write_Text);
           System.out.println(write_Text);
            }
            catch(Exception e){
                e.printStackTrace();
            }
    }
    
}


//package Declaration
package Panjala_Project1;

import java.util.ArrayList;

/**
 *
 * @author Uday Sagar Panjala - U00771807
 */
public class MainProgram {

    public static int k_features = 15;
    public static String file_name = "p1data.txt";

    public static void main(String[] a) {
        try {

            /// RUNNING 1ST TASK            
            ArrayList<Integer> list_of_genes_weka = new ArrayList<Integer>(k_features);
            WekaMethodTopKGenes obj = new WekaMethodTopKGenes();
            list_of_genes_weka = obj.WekaMethodTopKGenes1(k_features);

            /// RUNNING 2ND TASK            
            ArrayList<Integer> list_of_genes_entropy = new ArrayList<Integer>(k_features);
            EntropyBasedTopKGenes obj1 = new EntropyBasedTopKGenes();
            list_of_genes_entropy = obj1.EntropyBasedTopKGenes1(k_features);

            /// RUNNING 3RD TASK
            new CorrelationCoeff(file_name, k_features, list_of_genes_weka, list_of_genes_entropy);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

}

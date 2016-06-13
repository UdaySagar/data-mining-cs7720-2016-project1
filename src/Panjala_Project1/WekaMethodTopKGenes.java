//package Declaration
package Panjala_Project1;

import java.io.BufferedWriter;
import java.io.FileWriter;

import weka.core.Instances;
import weka.attributeSelection.CfsSubsetEval;
import weka.attributeSelection.BestFirst;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.supervised.attribute.AttributeSelection;
import java.util.*;
import java.util.List;
import java.util.Arrays;

public class WekaMethodTopKGenes {

    public ArrayList<Integer> WekaMethodTopKGenes1(int number_to_select) throws Exception {

        String dsfileName;//the data set file name
        Instances ini_instances = null; //the initial instances
        Instances slt_instances = null;   //the instances after attributes selection
        BufferedWriter out = null;
        int insnum;
        ArrayList<Integer> list_of_genes = new ArrayList<Integer>(number_to_select);

        dsfileName = System.getProperty("user.dir") + "/src/Panjala_Project1/cancerDataSet.arff";

        DataSource frData = new DataSource(dsfileName);
        ini_instances = frData.getDataSet();
        ini_instances.setClassIndex(ini_instances.numAttributes() - 1);

        AttributeSelection filter = new AttributeSelection();
        CfsSubsetEval cfseval = new CfsSubsetEval();
        BestFirst search = new BestFirst();
        filter.setEvaluator(cfseval);
        filter.setSearch(search);
        filter.setInputFormat(ini_instances);

        slt_instances = Filter.useFilter(ini_instances, filter);

        out = new BufferedWriter(new FileWriter(System.getProperty("user.dir") + "/src/Panjala_Project1/Output_WekaMethod.txt"));
   
        //output the data about the selected data set
        int k = 0;
        out.write("################################################################\n");
        out.write("Top K Genes selected using Weka Attribute Selection Method" + "\n");
        out.write("################################################################\n");
        System.out.println("###############################################################");
        System.out.println("Top K Genes selected using Weka Attribute Selection Method" + "\n");
        System.out.println("###############################################################");
        for (int j = 0; j < number_to_select; j++) {
            k = j + 1;
            out.write("Line no." + k + "--> Gene no. " + slt_instances.attribute(j).name() + "\n");
            System.out.println("Line no." + k + "--> Gene no. " + slt_instances.attribute(j).name() + "\n");
            list_of_genes.add(Integer.parseInt(slt_instances.attribute(j).name()));
        }
        out.close();
        return list_of_genes;
    }

}

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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Uday Sagar Panjala - U00771807
 */
public class EntropyBasedTopKGenes {

    HashMap<Integer, Double> hash = new HashMap<>();
    HashMap<Integer, Double> hash_for_geneInfogain = new HashMap<>();
    HashMap<Integer, Double> sorted_hash_for_geneInfogain = new HashMap<>();
    ArrayList lines = new ArrayList();
    static String[] fields;
    String[][] strings;
    Double[] arrayDouble;
    int left_final_count1;
    int left_final_count2;
    int left_final_count3;
    int right_final_count1;
    int right_final_count2;
    int right_final_count3;
    Double count_for_class1;
    Double count_for_class2;
    Double counts_split1;
    Double counts_split2;
    Double split1_pcount;
    Double split1_ncount;
    Double split2_pcount;
    Double split2_ncount;

    public ArrayList<Integer> EntropyBasedTopKGenes1(int k) throws FileNotFoundException, IOException {
        BufferedWriter out_entropy = new BufferedWriter(new FileWriter(System.getProperty("user.dir") + "/src/Panjala_Project1/Output_EntropyBased.txt"));
        ArrayList<Integer> list_of_genes1 = new ArrayList<Integer>(k);
        out_entropy.write("################################################################\n");
        out_entropy.write("Top K Genes selected using Entropy Based Method" + "\n");
        out_entropy.write("################################################################\n");
        System.out.println("###############################################################");
        System.out.println("Top K Genes selected using Entropy Based Method" + "\n");
        System.out.println("###############################################################");
        BufferedReader br = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("p1data.txt")));
        for (String line = br.readLine(); line != null; line = br.readLine()) {

            fields = line.split(", ");

            lines.add(fields);
        }
        strings = (String[][]) lines.toArray(new String[lines.size()][]);
        int i, j = 0;
        for (i = 0; i < fields.length; i++) {
            for (j = 0; j < strings.length; j++) {
                if (strings[j][i].contains(",")) {

                    String data_class[] = strings[j][i].split(",");
                    hash.put(j, Double.parseDouble(data_class[0]));

                } else {

                    hash.put(j, Double.parseDouble((strings[j][i])));

                }

            }
///// SENDING SINGLE ////////////////

            hash_for_geneInfogain.put(i, arrngeTheGenes(i));
        }
        sorted_hash_for_geneInfogain = sortByValues(hash_for_geneInfogain);

        int p = 0;
        int q = 0;
        for (Map.Entry<Integer, Double> entry : sorted_hash_for_geneInfogain.entrySet()) {
            p++;
            if (p < k + 1) {
                System.out.println("Line no. " + ++q + "--> Gene no. " + entry.getKey() + " -- Informatoin Gain: " + entry.getValue());
                out_entropy.write("Line no. " + ++q + "--> Gene no. " + entry.getKey() + " -- Informatoin Gain: " + entry.getValue() + "\n");
                list_of_genes1.add(entry.getKey());
            } else {
                break;
            }
        }

        out_entropy.close();
        return list_of_genes1;

    }

    public double arrngeTheGenes(int u) {
        arrayDouble = new Double[hash.size()];

        for (int first_array_index = 0; first_array_index < hash.size(); first_array_index++) {
            arrayDouble[first_array_index] = hash.get(first_array_index);

            for (int second_array_index = 0; second_array_index < hash.size(); second_array_index++) {
                if (arrayDouble[second_array_index] == null) {
                    arrayDouble[second_array_index] = 0.0;
                }
                if (arrayDouble[first_array_index] < arrayDouble[second_array_index]) {

                    Double value_to_be_swappped = arrayDouble[first_array_index];
                    arrayDouble[first_array_index] = arrayDouble[second_array_index];
                    arrayDouble[second_array_index] = value_to_be_swappped;

                }

            }
        }

        return spltTheGenes(u);

    }

    public double spltTheGenes(int g_num) {
        int no_of_intervals = 2;
        Double first_split = (arrayDouble[(arrayDouble.length - 1)] + arrayDouble[0]) / no_of_intervals;
        Double second_split_left = (first_split + arrayDouble[0]) / no_of_intervals;
        Double second_split_right = (first_split + arrayDouble[(arrayDouble.length - 1)]) / no_of_intervals;
        int count1_left = 0;
        int count2_left = 0;
        int count3_left = 0;
        int count1_right = 0;
        int count2_right = 0;
        int count3_right = 0;

////////////// LEFT SPLIT ////////////////////////////////         
        for (int z = 0; z < arrayDouble.length; z++) {

            if (arrayDouble[z] >= arrayDouble[0] && arrayDouble[z] <= second_split_left) {
                count1_left++;
                left_final_count1 = count1_left;

            }

            if (arrayDouble[z] > second_split_left && arrayDouble[z] <= first_split) {
                count2_left++;
                left_final_count2 = count2_left;

            }

            if (arrayDouble[z] > first_split && arrayDouble[z] <= arrayDouble[(arrayDouble.length - 1)]) {
                count3_left++;
                left_final_count3 = count3_left;

            }
        }

////////////// RIGHT SPLIT ////////////////////////////////         
        for (int z = 0; z < arrayDouble.length; z++) {

            if (arrayDouble[z] >= arrayDouble[0] && arrayDouble[z] <= second_split_left) {
                count1_right++;
                right_final_count1 = count1_right;

            }

            if (arrayDouble[z] > second_split_left && arrayDouble[z] <= first_split) {
                count2_right++;
                right_final_count2 = count2_right;

            }

            if (arrayDouble[z] > first_split && arrayDouble[z] <= arrayDouble[(arrayDouble.length - 1)]) {
                count3_right++;
                right_final_count3 = count2_right;

                String desc = "c, ";

            }
        }

        double info_gain_left = ComputeInfoGain(second_split_left, first_split);
        double info_gain_right = ComputeInfoGain(first_split, second_split_right);
        double info_gain;
        if (info_gain_left > info_gain_right) {
            info_gain = info_gain_left;
        } else {
            info_gain = info_gain_right;
        }

        return info_gain;

    }

    public double ComputeInfoGain(Double v1, Double v2) {
        Double gain = 0.0;
        Double count_for_leftbin = 0.0;
        Double count_for_rightbin = 0.0;
        Double p = 0.0;
        Double n = 0.0;
        Double p_for_firstbin = 0.0;
        Double n_for_firstbin = 0.0;
        Double p_for_secondbin = 0.0;
        Double n_for_secondbin = 0.0;
        Double p_for_thirdbin = 0.0;
        Double n_for_thirdbin = 0.0;
        Double mainEntrpy = 0.0;
        Double entropy_for_firstbin;
        Double entropy_for_secondbin;
        Double entropy_for_thirdbin;
        Double Info_Split;

        int i, j = 0;
        for (i = 0; i < fields.length; i++) {
            for (j = 0; j < strings.length; j++) {

                if (strings[j][i].contains(",")) {

                    String data_class[] = strings[j][i].split(",");

                    if (data_class[1].equalsIgnoreCase("positive")) {
                        String classVariable = data_class[1] + "\n";
                        Double gene = Double.parseDouble(strings[j][0]);
                        count_for_leftbin++;
                        count_for_class1 = count_for_leftbin;

                    } else if (data_class[1].equalsIgnoreCase("negative")) {
                        String classVariable = data_class[1] + "\n";
                        Double gene = Double.parseDouble(strings[j][0]);
                        count_for_rightbin++;
                        count_for_class2 = count_for_rightbin;

                    }

                }

            }

        }
        p = count_for_class1 / fields.length;
        n = count_for_class2 / fields.length;
        mainEntrpy = (-1) * ((p * Math.log10(p) / Math.log10(2.0)) + ((n)
                * Math.log10(n) / Math.log10(2.0)));

        Double split1counts1 = 0.0;
        Double p_counts_split1 = 0.0;
        Double n_counts_split1 = 0.0;
        Double split1counts2 = 0.0;
        Double p_counts_split2 = 0.0;
        Double n_counts_split2 = 0.0;

        Double split1counts3 = 0.0;
        Double v1_Counts3 = 0.0;
        Double pcounts3 = 0.0;
        Double p_count_s3 = 0.0;
        Double ncounts3 = 0.0;
        Double n_count_s3 = 0.0;

        for (int x = 0; x < arrayDouble.length; x++) {
            if (arrayDouble[x] <= v1) {
                split1counts1++;
                counts_split1 = split1counts1;
                if (strings[x][(fields.length - 1)].contains(",")) {

                    String data_class[] = strings[x][(fields.length - 1)].split(",");
                    if (data_class[1].equalsIgnoreCase("positive")) {

                        p_counts_split1++;
                        split1_pcount = p_counts_split1;

                    } else if (data_class[1].equalsIgnoreCase("negative")) {

                        n_counts_split1++;
                        split1_ncount = n_counts_split1;

                    }

                }
            } else if ((arrayDouble[x] > v1) && (arrayDouble[x] <= v2)) {
                split1counts2++;
                counts_split2 = split1counts2;
                if (strings[x][(fields.length - 1)].contains(",")) {

                    String data_class[] = strings[x][(fields.length - 1)].split(",");
                    if (data_class[1].equalsIgnoreCase("positive")) {

                        p_counts_split2++;
                        split2_pcount = p_counts_split2;

                    } else if (data_class[1].equalsIgnoreCase("negative")) {

                        n_counts_split2++;
                        split2_ncount = n_counts_split2;

                    }

                }
            } else {
                split1counts3++;
                v1_Counts3 = split1counts3;
                if (strings[x][(fields.length - 1)].contains(",")) {

                    String data_class[] = strings[x][(fields.length - 1)].split(",");

                    if (data_class[1].equalsIgnoreCase("positive")) {

                        pcounts3++;
                        p_count_s3 = pcounts3;

                    } else if (data_class[1].equalsIgnoreCase("negative")) {

                        ncounts3++;
                        n_count_s3 = ncounts3;

                    }

                }
            }
        }
        p_for_firstbin = split1_pcount / counts_split1;
        n_for_firstbin = split1_ncount / counts_split1;
        if (p_for_firstbin == 0) {
            entropy_for_firstbin = (-1) * (((n_for_firstbin) * Math.log10(n_for_firstbin) / Math.log10(2.0)));
        } else if (n_for_firstbin == 0) {
            entropy_for_firstbin = (-1) * ((p_for_firstbin * Math.log10(p_for_firstbin) / Math.log10(2.0)));
        } else {
            entropy_for_firstbin = (-1) * ((p_for_firstbin * Math.log10(p_for_firstbin) / Math.log10(2.0)) + ((n_for_firstbin)
                    * Math.log10(n_for_firstbin) / Math.log10(2.0)));
        }

        p_for_secondbin = split2_pcount / counts_split2;
        n_for_secondbin = split2_ncount / counts_split2;
        if (p_for_secondbin == 0) {
            entropy_for_secondbin = (-1) * (((n_for_secondbin) * Math.log10(n_for_secondbin) / Math.log10(2.0)));
        } else if (n_for_secondbin == 0) {
            entropy_for_secondbin = (-1) * ((p_for_secondbin * Math.log10(p_for_secondbin) / Math.log10(2.0)));
        } else {
            entropy_for_secondbin = (-1) * ((p_for_secondbin * Math.log10(p_for_secondbin) / Math.log10(2.0)) + ((n_for_secondbin)
                    * Math.log10(n_for_secondbin) / Math.log10(2.0)));
        }

        p_for_thirdbin = p_count_s3 / v1_Counts3;
        n_for_thirdbin = n_count_s3 / v1_Counts3;
        if (p_for_thirdbin == 0) {
            entropy_for_thirdbin = (-1) * (((n_for_thirdbin) * Math.log10(n_for_thirdbin) / Math.log10(2.0)));
        } else if (n_for_thirdbin == 0) {
            entropy_for_thirdbin = (-1) * ((p_for_thirdbin * Math.log10(p_for_thirdbin) / Math.log10(2.0)));
        } else {
            entropy_for_thirdbin = (-1) * ((p_for_thirdbin * Math.log10(p_for_thirdbin) / Math.log10(2.0)) + ((n_for_thirdbin)
                    * Math.log10(n_for_thirdbin) / Math.log10(2.0)));
        }

        Info_Split = (counts_split1 / fields.length) * (entropy_for_firstbin) + (split1counts2 / fields.length) * (entropy_for_secondbin) + (split1counts3 / fields.length) * (entropy_for_thirdbin);

        gain = mainEntrpy - Info_Split;

        return gain;
    }

    private static HashMap sortByValues(HashMap map) {
        List list = new LinkedList(map.entrySet());

        Collections.sort(list, new Comparator() {
            public int compare(Object o1, Object o2) {
                return ((Comparable) ((Map.Entry) (o2)).getValue())
                        .compareTo(((Map.Entry) (o1)).getValue());
            }
        });
        HashMap sortedHashMap = new LinkedHashMap();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Map.Entry entry = (Map.Entry) it.next();
            sortedHashMap.put(entry.getKey(), entry.getValue());
        }
        return sortedHashMap;
    }
}

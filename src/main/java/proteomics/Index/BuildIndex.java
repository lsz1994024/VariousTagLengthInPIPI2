/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Index;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMIndex;
import proteomics.FM.UnicodeReader;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.*;
import proteomics.Types.PepInfo;
//import org.apache.commons.langs.ArrayUtils;


public class BuildIndex {
    public final MassTool massTool;
    private Map<Character, Double> fixModMap = new HashMap<>(25, 1);
    public double minPeptideMass = 9999;
    public double maxPeptideMass = 0;
    public InferSegment inferSegment;
    public TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();
    public Map<String, PepInfo> pepInfoMap;
    private final String labelling;
    public final DbTool dbTool; // this one doesn't contain contaminant proteins.
    private InferPTM inferPTM;
    public Map<String, Integer> protLengthMap = new HashMap<>();
    public FMIndex fmIndexFull;
    public int[] dotPosArrFull;
    public Map<Integer, String> posProtMapFull = new HashMap<>();
    public Map<String, String> protSeqMap;
    public Map<String, Set<Pair<String,Integer>>> tagProtPosMap = new HashMap<>();
    public BuildIndex(Map<String, String> parameterMap) throws Exception {
        boolean addContaminant = parameterMap.get("add_contaminant").contentEquals("1");
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = parameterMap.get("15N").trim().contentEquals("1") ? "N15" : "N14";

        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('X', 0.0);  //X for any amino acid, only used for NC term mod on any site, specified by user
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Double.valueOf(parameterMap.get("c")));

        // read protein database
        dbTool = new DbTool(dbPath, parameterMap.get("database_type"));
        DbTool contaminantsDb = null;
        if (addContaminant) {
            contaminantsDb = new DbTool(null, "contaminants");
            protSeqMap = contaminantsDb.getProtSeqMap();
            protSeqMap.putAll(dbTool.getProtSeqMap()); // using the target sequence to replace contaminant sequence if there is conflict.
        } else {
            protSeqMap = dbTool.getProtSeqMap();
        }
//        DECOY_Q29443
        // define a new MassTool object
        massTool = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site_1").trim(), parameterMap.get("protection_site_1").trim(), parameterMap.get("is_from_C_term_1").trim().contentEquals("1"), parameterMap.getOrDefault("cleavage_site_2", null), parameterMap.getOrDefault("protection_site_2", null), parameterMap.containsKey("is_from_C_term_2") ? parameterMap.get("is_from_C_term_2").trim().contentEquals("1") : null, ms2Tolerance, oneMinusBinOffset, labelling);

        inferPTM = new InferPTM(massTool, fixModMap, parameterMap);

        // build database
        inferSegment = new InferSegment(massTool, parameterMap, fixModMap);

        BufferedWriter writerProt = new BufferedWriter(new FileWriter("catProt.txt"));
        int dotPos = 0;
        int dotNum = 0;
        dotPosArrFull = new int[protSeqMap.keySet().size()];

        for (String protId : protSeqMap.keySet()) {
            dotPosArrFull[dotNum] = dotPos;
            posProtMapFull.put(dotNum, protId);
            String protSeq = protSeqMap.get(protId).replace('I', 'L');
            protSeqMap.put(protId, protSeq);
            writerProt.write("." + protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length()+1;
            int numOfTags = inferSegment.getLongTagNumForProt(protSeq);
            protLengthMap.put(protId, numOfTags);
        }
        writerProt.close();
        char[] text = loadFile("catProt.txt", true);
        fmIndexFull = new FMIndex(text);
        System.out.println("Finish build FM index");


    }

    public static char[] loadFile(String file, boolean appendTerminalCharacter) throws IOException{
        // read text file, auto recognize bom marker or use
        // system default if markers not found.
        BufferedReader reader = null;
        CharArrayWriter writer = null;
        UnicodeReader r = new UnicodeReader(new FileInputStream(file), null);

        char[] buffer = new char[16 * 1024];   // 16k buffer
        int read;
        try {
            reader = new BufferedReader(r);
            writer = new CharArrayWriter();
            while( (read = reader.read(buffer)) != -1) {
                writer.write(buffer, 0, read);
            }
            if (appendTerminalCharacter) {
                writer.append('\0');
            }
            writer.flush();
            return writer.toCharArray();
        } catch (IOException ex) {
            throw ex;
        } finally {
            try {
                writer.close(); reader.close(); r.close();
            } catch (Exception ex) { }
        }
    }
    public DbTool getDbTool() {
        return dbTool;
    }

    public MassTool returnMassTool() {
        return massTool;
    }

    public double getMinPeptideMass() {
        return minPeptideMass;
    }

    public double getMaxPeptideMass() {
        return maxPeptideMass;
    }

    public Map<Character, Double> returnFixModMap() {
        return fixModMap;
    }

    public InferSegment getInferSegment() {
        return inferSegment;
    }

    public InferPTM getInferPTM() {
        return inferPTM;
    }

    public TreeMap<Double, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    public Map<String, PepInfo> getPepInfoMap() {
        return pepInfoMap;
    }

    public String getLabelling() {
        return labelling;
    }


}

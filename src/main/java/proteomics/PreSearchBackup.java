///*
// * Copyright 2016-2019 The Hong Kong University of Science and Technology
// *
// * Licensed under the Apache License, Version 2.0 (the "License");
// * you may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *     http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */
//
//package proteomics;
//
//import ProteomicsLibrary.MassTool;
//import ProteomicsLibrary.SpecProcessor;
//import org.apache.commons.math3.util.Pair;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;
//import proteomics.Index.BuildIndex;
//import proteomics.Search.Search;
//import proteomics.Segment.InferSegment;
//import proteomics.Spectrum.DatasetReader;
//import proteomics.Types.ExpTag;
//import proteomics.Types.Peptide;
//import proteomics.Types.PeptideInfo;
//import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
//
//import java.util.*;
//import java.util.concurrent.Callable;
//import java.util.concurrent.locks.ReentrantLock;
//
//import static proteomics.PIPI.lszDebugScanNum;
//
//public class PreSearchBackup implements Callable<PreSearchBackup.Entry> {
//    private static final Logger logger = LoggerFactory.getLogger(PreSearchBackup.class);
//
//    private static final int candisNum = 20;
//    private final BuildIndex buildIndex;
//    private final MassTool massTool;
//    private final double ms1Tolerance;
//    private final double leftInverseMs1Tolerance;
//    private final double rightInverseMs1Tolerance;
//    private final int ms1ToleranceUnit;
//    private final double minPtmMass;
//    private final double maxPtmMass;
//    private final int localMaxMs2Charge;
//    private final JMzReader spectraParser;
//    private final double minClear;
//    private final double maxClear;
//    private final ReentrantLock lock;
//    private final String scanName;
//    private final int precursorCharge;
//    private final double precursorMass;
//    private final SpecProcessor specProcessor;
//    private final int scanNum;
//    private String truth;
//
//
//    public PreSearchBackup(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
//            , int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
//            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
//            , SpecProcessor specProcessor, String truth) {
//
//        this.buildIndex = buildIndex;
//        this.massTool = massTool;
//        this.ms1Tolerance = ms1Tolerance;
//        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
//        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
//        this.ms1ToleranceUnit = ms1ToleranceUnit;
//        this.minPtmMass = minPtmMass;
//        this.maxPtmMass = maxPtmMass;
//        this.localMaxMs2Charge = localMaxMs2Charge;
//        this.spectraParser = spectraParser;
//        this.minClear = minClear;
//        this.maxClear = maxClear;
//        this.lock = lock;
//        this.scanName = scanName;
//        this.precursorCharge = precursorCharge;
//        this.precursorMass = precursorMass;
//        this.specProcessor = specProcessor;
//        this.scanNum = scanNum;
//        this.truth = truth;
//    }
//
//    @Override
//    public Entry call() throws Exception {
//        Map<Double, Double> rawPLMap;
//        try {
//            lock.lock();
//            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
//        } finally {
//            lock.unlock();
//        }
//        if (lszDebugScanNum.contains(this.scanNum)) {
//            int a = 1;
//        }
//        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);
//
//        if (plMap.isEmpty()) return null;
//        // Coding
//        InferSegment inferSegment = buildIndex.getInferSegment();
//        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
//
//        List<ExpTag> tag4List = inferSegment.getAllTag4(precursorMass, finalPlMap, scanNum);
//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, 4);
//
//        Collections.sort(tag4List, Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
//
//        String truth = "DNVFENNRLAFEVAEK";
//        if (lszDebugScanNum.contains(scanNum)){
//            for (ExpTag tagInfo : tag4List){
//                String seq = tagInfo.getFreeAaString();
//                String revSeq = (new StringBuilder(seq)).reverse().toString();
//                if (truth.contains(seq) || truth.contains(revSeq)) {
//                    System.out.println(seq +"  1");
//                }
//            }
//        }
//        Map<String, Set<Pair<String, Integer>>> tagProtPosMap = buildIndex.tagProtPosMap;
//        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);
//
//        if (tag4List.isEmpty())  return null;
//        Entry entry = new Entry();
//
//
//        double totalMass = precursorMass + 2 * MassTool.PROTON;
//        List<ExpTag> compTag4List = new ArrayList<>(2* Math.min(75, tag4List.size())); // only take top 100 candidates.
//        compTag4List.addAll(tag4List.subList(0, Math.min(75, tag4List.size())));
//        for (ExpTag tag4 : tag4List.subList(0, Math.min(75, tag4List.size()))) {
//            compTag4List.add(tag4.revTag(totalMass));
//        }   // duplicate reverse tags for every tag, then dont need to care the both direction of a tag because already contained in this compTag4List
//
////        logger.info(scanNum+",starts,"+compTag4List.size());
//        double pcMassL = precursorMass - 250;
//        double pcMassR = precursorMass + 250;
//        String tag;
//        for (ExpTag tagInfo : compTag4List) {
//            tag = tagInfo.getFreeAaString();
//
//            if (!tagProtPosMap.containsKey(tagInfo.getFreeAaString())) {
//                int a = 1;
//                continue;
//            }
//            Set<Pair<String, Integer>> protPosPairs = tagProtPosMap.get(tagInfo.getFreeAaString());
//            for (Pair<String, Integer> protPos : protPosPairs){
//                String protId = protPos.getFirst();
//                double tagCMass = tagInfo.getTailLocation();
//                int pos = protPos.getSecond();
//                String protSeq = buildIndex.protSeqMap.get(protId);
//                Set<Integer> cPosSet = new HashSet<>();
//
//                if (isKR(protSeq.charAt(pos+4)) && tagCMass <= pcMassR && tagCMass  >= pcMassL) {
//                    cPosSet.add(pos+4);
//                }
//                for (int i = pos+tag.length(); i < protSeq.length(); i++) {  //
//                    if (isX(protSeq.charAt(i))) break;
//
//                    tagCMass += massTool.getMassTable().get(protSeq.charAt(i));
//                    if (tagCMass < pcMassL) continue;
//                    if (tagCMass > pcMassR) break;
//                    if (isKR(protSeq.charAt(i))) {
//                        cPosSet.add(i);
//                    }
//                }
//                for (int cPos : cPosSet) {
//                    double deltaMass = precursorMass - massTool.calResidueMass(protSeq.substring(pos,cPos+1)) - massTool.H2O;
//                    char rightFlank;
//                    if (cPos == protSeq.length()-1) {
//                        rightFlank = '-';
//                    } else {
//                        rightFlank = protSeq.charAt(cPos+1);
//                    }
//
////                    int numMissCleave = 0;
//                    for (int nPos = pos-1; nPos > 0; nPos--) {
//                        if (isX(protSeq.charAt(nPos))) break;
//                        deltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
//                        if (deltaMass > 250) continue;
//                        if (deltaMass < -250) break;
//                        String pepSeq = protSeq.substring(nPos, cPos+1);
//
//                        char leftFlank;
//                        if (nPos == 0 || (nPos == 1 && protSeq.charAt(0) == 'M')){
//                            leftFlank = '-';
//                        } else {
//                            leftFlank = protSeq.charAt(nPos-1);
//                        }
//
//                        if (peptideInfoMap.containsKey(pepSeq)) {
//                            PeptideInfo pepInfo = peptideInfoMap.get(pepSeq);
//                            if (pepInfo.leftFlank != '-' && pepInfo.rightFlank != '-') {
//                                if (rightFlank == '-' || leftFlank == '-') {
//                                    pepInfo.leftFlank = leftFlank;
//                                    pepInfo.rightFlank = rightFlank;
//                                }
//                            }
//                            pepInfo.protIdSet.add(protId);
//                            if (!protId.startsWith("DECOY_")) {
//                                pepInfo.isTarget = true;
//                            }
//                        } else {
//                            PeptideInfo pepInfo = new PeptideInfo(pepSeq, !protId.startsWith("DECOY_"), leftFlank, rightFlank);
//                            pepInfo.protIdSet.add(protId);
//                            peptideInfoMap.put(pepSeq, pepInfo);
//                        }
//                    }
//                }
//            }
//        }
//        entry.scanName = this.scanName;
////        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, 4);
//
//        Map<String, Double> scanTagStrMap = inferSegment.getTagStrMap(allLongTagList);
//        Search search = new Search(entry, scanNum, buildIndex.inferSegment, precursorMass, scanTagStrMap, massTool, localMaxMs2Charge, peptideInfoMap);
////        entry.peptideInfoMap = peptideInfoMap;
//        if (lszDebugScanNum.contains(scanNum)){
//            if (peptideInfoMap.containsKey(truth)) {
//                int a = 1;
//            }
//        }
//        for (Peptide pep : entry.ptmOnlyList ) {
//            entry.peptideInfoMap.put(pep.getPTMFreePeptide(), peptideInfoMap.get(pep.getPTMFreePeptide()).clone());
//        }
//        for (Peptide pep : entry.ptmFreeList ) {
//            entry.peptideInfoMap.put(pep.getPTMFreePeptide(), peptideInfoMap.get(pep.getPTMFreePeptide()).clone());
//        }
//        int c = 1;
//        return entry;
//    }
//
//    private boolean isKR(char aa){
//        return aa == 'K' || aa == 'R';
//    }
//    private boolean isX(char aa){
//        return aa == 'X';
//    }
//    public class Entry {
//
//        public Map<String, PeptideInfo> peptideInfoMap = new HashMap<>();
//        public double precursorMass = PreSearchBackup.this.precursorMass;
//        public int precursorCharge = PreSearchBackup.this.precursorCharge;
//        public List<Peptide> ptmOnlyList = new ArrayList<>();
//        public List<Peptide> ptmFreeList = new ArrayList<>();
//
//        public String scanName;
//        Entry() {
//        }
//    }
//}

package proteomics.FM;

import java.io.Serializable;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;

public class FMIndex implements Serializable {

    public static long space = 0;

    // text
    public char[] S;
    // text length
    public int m;
    // suffix array
    public Integer[] SA;
    // Burrows-Wheeler Transform of text
    public ArrayList<Character> BWT;

    // C.get(c) C(c) returns the number of occurrences in S of characters alphabetically smaller than c
    public HashMap<Character, Integer> C;
    public ArrayList<Character> sortedAlphabet;
    public int sigma;

    public WaveletTreeNode waveletTreeRoot;

//    @Override
//    public boolean equals(Object o) {
//        FMIndex c = (FMIndex) o;
//        boolean res = true;
//        if (space != c.space) return false;
//        for (char cc : S){
//
//        }
//        return
//    }


    public FMIndex(char[] text){
        S = text ;
        m = text.length;

        createSuffixArray();

        createBWT();

        // create wavelet tree of BWT
        createWaveletTree();

    }

    public SearchInterval fmSearchFuzzy(char[] P){
        int i = P.length - 1;
        int sp = 0;
        int ep = m - 1;
//        int absTagPos = fmIndex.SA[ii];
        int lastSp = sp;
        int lastEp = ep;
        while (sp <= ep && i >= 0) {
            if (!C.containsKey(P[i])) {
                return null;
            }
            lastSp = sp;
            lastEp = ep;
            sp = C.get(P[i]) + Occ(sp - 1, P[i]);
            ep = C.get(P[i]) + Occ(ep, P[i]) - 1;
            i --;
        }

        if (ep < sp) {
            SearchInterval unSettledSI = new SearchInterval(lastSp,lastEp);
            unSettledSI.settled = false;
            unSettledSI.matchedPos = i + 2;
            return unSettledSI;
        } else {
            return new SearchInterval(sp, ep);
        }
    }

    public SearchInterval fmSearch(char[] P){
        int i = P.length - 1;
        int sp = 0;
        int ep = m - 1;
//        int absTagPos = fmIndex.SA[ii];
        while (sp <= ep && i >= 0) {
            if (!C.containsKey(P[i])) {
                return null;
            }
            sp = C.get(P[i]) + Occ(sp - 1, P[i]);
            ep = C.get(P[i]) + Occ(ep, P[i]) - 1;
            i --;
        }

        if (ep < sp) {
            return null;
        } else {
            return new SearchInterval(sp, ep);
        }
    }


    // create suffix array using Manber and Myers algorithm
    private void createSuffixArray() {

        C = new HashMap<Character, Integer>();
        sortedAlphabet = new ArrayList<Character>();

        Integer[] SAH = new Integer[m];
        for ( int i = 0; i < m; i++ ){
            SAH[i] = i;
        }
        Integer[] SA2H = new Integer[m];
        int[] bucketID = new int[m];
        int[] newBucketID = new int[m];
        int[] bucketStartIndex = new int[m];
        int[] bucketCounter = new int[m];
        int[] bucketTracking = new int[m];

        // initialize SA1
        IndexSorter is = new IndexSorter(S, SAH);
        is.sort();

        int numberOfBuckets = 1;
        bucketTracking[0] = -1;
        sortedAlphabet.add(S[SAH[0]]);
        // initialize buckets data structures
        for (int i = 1; i < m; i++) {
            if (S[SAH[i]] == S[SAH[i - 1]]) {
                // still in the same bucket
                bucketID[SAH[i]] = bucketID[SAH[i - 1]];
            } else {
                // start of new bucket
                C.put(S[SAH[i]], i);
                sortedAlphabet.add(S[SAH[i]]);
                bucketID[SAH[i]] = bucketID[SAH[i - 1]] + 1;
                bucketStartIndex[bucketID[SAH[i]]] = i;
                numberOfBuckets++;
            }
            bucketTracking[i] = -1;
        }

        sigma = C.size();

        try {
            Field tableField = HashMap.class.getDeclaredField("table");
            tableField.setAccessible(true);
            Object[] table;
            table = (Object[]) tableField.get(C);

            FMIndex.space += 32 * C.size() +  4 * (table == null ? 0 : table.length);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int j, p;
        int H = 1;

        while (numberOfBuckets < m) {
            for (int i = 0; i < m; i++) {
                j = SAH[i];

                // Aj is in a singelton bucket if start index of the buket is the last index,
                // or if start index of the bucket only differs by one from starting index of the next bucket
                if (bucketStartIndex[bucketID[j]] == m - 1 || bucketStartIndex[bucketID[j]] == bucketStartIndex[bucketID[j] + 1] - 1) {
                    // update only if it's the first time
                    if (bucketTracking[j] == -1) {
                        SA2H[i] = j;
                        bucketTracking[j] = bucketID[j];
                    }
                }

                if (j - H >= 0) {
                    // update only if it's the first time
                    if (bucketTracking[j - H] == -1) {
                        p = bucketStartIndex[bucketID[j - H]] + bucketCounter[bucketID[j - H]];
                        SA2H[p] = j - H;
                        bucketCounter[bucketID[j - H]] ++;
                        bucketTracking[j - H] = bucketID[j];
                    }
                }
            }
            // update new buckets
            bucketTracking[SA2H[0]] = -1; // tracking of null character
            numberOfBuckets = 1;
            for (int i = 1; i < m; i++) {
                if ((bucketID[SAH[i]] == bucketID[SAH[i - 1]]) && bucketTracking[SA2H[i]] == bucketTracking[SA2H[i - 1]]) {
                    // still in the same bucket
                    newBucketID[SA2H[i]] = newBucketID[SA2H[i - 1]];
                } else {
                    // start of new bucket
                    newBucketID[SA2H[i]] = newBucketID[SA2H[i - 1]] + 1;
                    bucketStartIndex[newBucketID[SA2H[i]]] = i;
                    numberOfBuckets++;
                }
            }

            // reset and copy to prepare for next loop
            for (int i = 0; i < m; i++) {
                // copy newBucketID to bucketID
                bucketID[SA2H[i]] = newBucketID[SA2H[i]];
                // reset bucket tracking
                bucketTracking[SA2H[i]] = -1;
                // copy new SA2H to SAH
                SAH[i] = SA2H[i];
                // reset bucket counter
                bucketCounter[i] = 0;
            }

            H = H * 2;
        }

        SA = SAH;

        FMIndex.space += 16 * SA.length;
    }

    // get Burrows-Wheeler Transform from suffix array
    private void createBWT() {
        BWT = new ArrayList<Character>(m);

        for (int i = 0; i < m; i++) {
            if (SA[i] == 0) {
                BWT.add(S[m - 1]);
            } else {
                BWT.add(S[SA[i] - 1]);
            }
        }
    }

    private void createWaveletTree(){
        waveletTreeRoot = new WaveletTreeNode(BWT, sortedAlphabet);
    }

    private int Occ(int i, char c){
        return waveletTreeRoot.Occ(i + 1, c);
    }
}

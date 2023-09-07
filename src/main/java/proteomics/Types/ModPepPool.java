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

package proteomics.Types;

import java.util.Comparator;
import java.util.TreeSet;

public class ModPepPool {
    private int maxNumOfPtmPatterns = 10;

    public final String freeSeq;
    public Peptide bestPep = null;
    public TreeSet<Peptide> peptideTreeSet = new TreeSet<>(Comparator.reverseOrder());

    public ModPepPool(String freeSeq) {
        this.freeSeq = freeSeq;
    }

    public ModPepPool(String freeSeq, int maxNumOfPtmPatterns) {
        this.freeSeq = freeSeq;
        this.maxNumOfPtmPatterns = maxNumOfPtmPatterns;
    }
    public Peptide getTopPepPtn() {
        return peptideTreeSet.first();
    }

    public void push(Peptide peptide) {
        if (peptideTreeSet.size() < maxNumOfPtmPatterns) { //max restore 5 patterns for one peptide
            peptideTreeSet.add(peptide);
        } else if (peptideTreeSet.last().compareTo(peptide) < 0) {
            peptideTreeSet.pollLast();
            peptideTreeSet.add(peptide);
        }
    }

    public TreeSet<Peptide> getPeptideTreeSet() {
        return peptideTreeSet;
    }
}

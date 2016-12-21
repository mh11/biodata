/*
 * Copyright 2015 OpenCB
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

package org.opencb.biodata.models.variant;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.StringUtils;
import org.opencb.biodata.models.variant.avro.StructuralVariation;
import org.opencb.biodata.models.variant.avro.VariantAnnotation;
import org.opencb.biodata.models.variant.avro.VariantAvro;
import org.opencb.biodata.models.variant.avro.VariantType;

import java.io.Serializable;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

/**
 * @author Jacobo Coll;
 * @author Cristina Yenyxe Gonzalez Garcia &lt;cyenyxe@ebi.ac.uk&gt;
 */
@JsonIgnoreProperties({"impl", "ids", "sourceEntries", "studiesMap", "lengthReference", "lengthAlternate"})
public class Variant implements Serializable, Comparable<Variant> {

    public static final EnumSet<VariantType> SV_SUBTYPES = EnumSet.of(VariantType.INSERTION, VariantType.DELETION,
            VariantType.TRANSLOCATION, VariantType.INVERSION, VariantType.CNV);
    private final AtomicReference<VariantAvro> impl = new AtomicReference<>();
    private volatile Map<String, StudyEntry> studyEntries = null;

    public static final int SV_THRESHOLD = 50;
    public static final String CNVSTR = "<CN";

    public Variant() {
        impl.set(new VariantAvro(null, new LinkedList<>(), "", -1, -1, "", "", "+", null, 0, null, new HashMap<>(), new LinkedList<>(), null));
    }

    public Variant(VariantAvro avro) {
        Objects.requireNonNull(avro);
        impl.set(avro);
    }

    /**
     * Creates a variant parsing a string.
     *
     * Format : (chr):(start)[-(end)][:(ref)]:(alt)
     *
     * @param variantString Variant string
     * @throws IllegalArgumentException if the variant does not match with the pattern
     */
    public Variant(String variantString) {
        this();
        if (variantString != null && !variantString.isEmpty()) {
            String[] fields = variantString.split(":", -1);
            if (fields.length == 3) {
                setChromosome(fields[0]);
                setAlternate(checkEmptySequence(fields[2]));
                // Structural variant needs start-end coords
                if (fields[1].contains("-")) {
                    String[] coordinatesParts = fields[1].split("-");
                    setReference("N");
                    setStart(Integer.parseInt(coordinatesParts[0]));
                    setEnd(Integer.parseInt(coordinatesParts[1]));
                    setLength(inferLengthSV(getAlternate(), getStart(), getEnd()));
                // Short variant, no reference specified
                } else {
                    setStart(Integer.parseInt(fields[1]));
                    setReference("");
                    setLength(inferLengthShortVariant(getReference(), getAlternate()));
                    setEnd(getStart() + getLengthReference() - 1);
                }
            } else {
                if (fields.length == 4) {
                    setChromosome(fields[0]);
                    setAlternate(checkEmptySequence(fields[3]));
                    // Structural variant needs start-end coords
                    if (fields[1].contains("-")) {
                        String[] coordinatesParts = fields[1].split("-");
                        setReference(checkEmptySequence(fields[2]));
                        setStart(Integer.parseInt(coordinatesParts[0]));
                        setEnd(Integer.parseInt(coordinatesParts[1]));
                        setLength(inferLengthSV(getAlternate(), getStart(), getEnd()));
                    } else {
                        setStart(Integer.parseInt(fields[1]));
                        setReference(checkEmptySequence(fields[2]));
                        setLength(inferLengthShortVariant(getReference(), getAlternate()));
                        setEnd(getStart() + getLengthReference() - 1);
                    }
                } else {
                    throw new IllegalArgumentException("Variant " + variantString + " needs 3 or 4 fields separated by ':'. "
                            + "Format: \"(chr):(start)[-(end)][:(ref)]:(alt)\"");
                }
            }
        }
        resetType();

        if (VariantType.CNV.equals(getType())) {
            setSv(new StructuralVariation(this.getStart(), this.getStart(), this.getEnd(), this.getEnd(),
                    getCopyNumberFromStr(this.getAlternate())));

        }
    }

    public Variant(String chromosome, int position, String reference, String alternate) {
        this(chromosome, position, position, reference, alternate, "+");
        setEnd(getStart() + getLengthReference() - 1);
    }

    public Variant(String chromosome, int start, int end, String reference, String alternate) {
        this(chromosome, start, end, reference, alternate, "+");
    }

    public Variant(String chromosome, int start, int end, String reference, String alternate, String strand) {
        impl.set(new VariantAvro(
                null,
                new LinkedList<>(),
                "",
                start,
                end,
                checkEmptySequence(reference),
                checkEmptySequence(alternate),
                strand,
                null,
                0,
                null,
                new HashMap<>(),
                new LinkedList<>(),
                null));
        if (start > end && !(reference.equals("-") || reference.isEmpty())) {
            throw new IllegalArgumentException("End position must be greater than the start position for variant: "
                    + chromosome + ":" + start + "-" + end + ":" + reference + ":" + alternate);
        }

        this.setChromosome(chromosome);

        this.resetLength();
        this.resetType();

        if (VariantType.CNV.equals(getType())) {
            setSv(new StructuralVariation(this.getStart(), this.getStart(), this.getEnd(), this.getEnd(),
                    getCopyNumberFromStr(this.getAlternate())));

        }

//        this.resetHGVS();

//        this.annotation = new VariantAnnotation(this.chromosome, this.start, this.end, this.reference, this.alternate);
        studyEntries = new HashMap<>();
    }

    private Integer getCopyNumberFromStr(String cnvStr) {
        String copyNumberString = cnvStr.split(CNVSTR)[1].split(">")[0];
        if (StringUtils.isNumeric(copyNumberString)) {
            return Integer.valueOf(copyNumberString);
        } else {
            return null;
        }
    }

    private String checkEmptySequence(String sequence) {
        return (sequence != null && !sequence.equals("-")) ? sequence : "";
    }

    private void resetType() {
        setType(inferType(getReference(), getAlternate(), getLength()));
    }

    public static VariantType inferType(String reference, String alternate, Integer length) {
        if (Allele.wouldBeSymbolicAllele(alternate.getBytes()) || Allele.wouldBeSymbolicAllele(reference.getBytes())) {
            if (alternate.startsWith(CNVSTR)) {
                return VariantType.CNV;
            } else {
                return VariantType.SYMBOLIC;
            }
        } else if (alternate.equals(Allele.NO_CALL_STRING)) {
            return VariantType.NO_VARIATION;
        } else {
            if (reference.length() == alternate.length()) {
                if (length > 1) {
                    return VariantType.MNV;
                } else {
                    return VariantType.SNV;
                }
            } else {
                if (length <= SV_THRESHOLD) {
                /*
                * 3 possibilities for being an INDEL:
                * - The value of the ALT field is <DEL> or <INS>
                * - The REF allele is not . but the ALT is
                * - The REF allele is . but the ALT is not
                * - The REF field length is different than the ALT field length
                */
                    return VariantType.INDEL;
                } else {
                    return VariantType.SV;
                }
            }
        }
    }

    public void resetLength() {
        setLength(inferLength(getReference(), getAlternate(), getStart(), getEnd()));
    }

    public static int inferLength(String reference, String alternate, int start, int end) {
        final int length;
        if (reference == null) {
            length = inferLengthSV(alternate, start, end);
        } else {
            length = inferLengthShortVariant(reference, alternate);
        }
        return length;
    }

    private static int inferLengthShortVariant(String reference, String alternate) {
        final int length;
        if (alternate == null) {
            length = reference.length();
        } else {
            length = Math.max(reference.length(), alternate.length());
        }
        return length;
    }

    private static int inferLengthSV(String alternate, int start, int end) {
        int length;
        if (StringUtils.startsWith(alternate, CNVSTR)) {
            length = end - start + 1;
        } else {
            length = alternate == null ? 0 : alternate.length();
        }
        return length;
    }

//    public void resetHGVS() {
//        if (this.getType() == VariantType.SNV || this.getType() == VariantType.SNP) { // Generate HGVS code only for SNVs
//            List<String> hgvsCodes = new LinkedList<>();
//            hgvsCodes.add(getChromosome() + ":g." + getStart() + getReference() + ">" + getAlternate());
//            if (impl.getHgvs() == null) {
//                impl.setHgvs(new HashMap<>());
//            }
//            impl.getHgvs().put("genomic", hgvsCodes);
//        }
//    }

    public VariantAvro getImpl() {
        return impl.get();
    }

    public final void setChromosome(String chromosome) {
        if (chromosome == null || chromosome.length() == 0) {
            throw new IllegalArgumentException("Chromosome must not be empty");
        }
        // Replace "chr" references only at the beginning of the chromosome name
        // For instance, tomato has SL2.40ch00 and that should be kept that way
        if (chromosome.startsWith("ch")) {
            if (chromosome.startsWith("chrom")) {
                getImpl().setChromosome(chromosome.substring(5));
            } else if (chromosome.startsWith("chrm")) {
                getImpl().setChromosome(chromosome.substring(4));
            } else if (chromosome.startsWith("chr")) {
                getImpl().setChromosome(chromosome.substring(3));
            } else {
                // Only starts with ch
                getImpl().setChromosome(chromosome.substring(2));
            }
        } else {
            getImpl().setChromosome(chromosome);
        }
    }

    public final void setStart(Integer start) {
        if (start < 0) {
            throw new IllegalArgumentException("Start must be positive");
        }
        getImpl().setStart(start);
    }

    public final void setEnd(Integer end) {
        if (end < 0) {
            throw new IllegalArgumentException("End must be positive");
        }
        getImpl().setEnd(end);
    }

    public void setReference(String reference) {
        getImpl().setReference(reference);
//        resetLength();
    }

    public void setAlternate(String alternate) {
        getImpl().setAlternate(alternate);
//        resetLength();
    }

    public String getId() {
        return getImpl().getId();
    }

    public void setId(String id) {
        getImpl().setId(id);
    }

    public String getChromosome() {
        return getImpl().getChromosome();
    }

    public Integer getStart() {
        return getImpl().getStart();
    }

    public Integer getEnd() {
        return getImpl().getEnd();
    }

    public String getReference() {
        return getImpl().getReference();
    }

    public String getAlternate() {
        return getImpl().getAlternate();
    }

    public String getStrand() {
        return getImpl().getStrand();
    }

    public void setStrand(String strand) {
        getImpl().setStrand(strand);
    }

    public StructuralVariation getSv() {
        return getImpl().getSv();
    }

    public void setSv(StructuralVariation sv) {
        getImpl().setSv(sv);
    }

    @Deprecated
    public List<String> getIds() {
        if (StringUtils.isNotEmpty(getImpl().getId())) {
            if (getImpl().getNames() != null) {
                List<String> ids = new ArrayList<>(1 + getImpl().getNames().size());
                ids.add(getImpl().getId());
                ids.addAll(getImpl().getNames());
                return ids;
            } else {
                return Collections.singletonList(getImpl().getId());
            }
        } else {
            return getImpl().getNames();
        }
    }

    @Deprecated
    public void setIds(List<String> ids) {
        if (ids == null || ids.isEmpty()) {
            getImpl().setId(null);
            getImpl().setNames(Collections.emptyList());
        } else {
            getImpl().setId(ids.get(0));
            getImpl().setNames(ids.subList(1, ids.size()));
        }
    }

    public List<String> getNames() {
        return getImpl().getNames();
    }

    public void setNames(List<String> names) {
        getImpl().setNames(names);
    }

    public Integer getLength() {
        return getImpl().getLength();
    }

    public Integer getLengthReference() {
        if (EnumSet.of(VariantType.NO_VARIATION, VariantType.CNV, VariantType.SV, VariantType.SYMBOLIC).contains(getType())) {
            return getLength();
        } else {
            return getReference().length();
        }
    }

    public Integer getLengthAlternate() {
        return getAlternate().length();
    }

    public void setLength(Integer value) {
        getImpl().setLength(value);
    }

    public VariantType getType() {
        return getImpl().getType();
    }

    public void setType(VariantType value) {
        getImpl().setType(value);
    }

    public Map<String, List<String>> getHgvs() {
        return getImpl().getHgvs();
    }

    public void setHgvs(Map<String, List<String>> value) {
        getImpl().setHgvs(value);
    }

    public VariantAnnotation getAnnotation() {
        return getImpl().getAnnotation();
    }

    public void setAnnotation(VariantAnnotation value) {
        getImpl().setAnnotation(value);
    }

    public boolean addHgvs(String type, String value) {
        List<String> listByType = getHgvs().get(type);
        if (listByType == null) {
            listByType = new LinkedList<>();
        }
        if (!listByType.contains(value)) {
            return listByType.add(value);
        } else {
            return false; //Collection has not changed
        }
    }

    public List<StudyEntry> getStudies() {
        return getStudiesMap() == null ? null : Collections.unmodifiableList(new ArrayList<>(getStudiesMap().values()));
    }

    public void setStudies(List<StudyEntry> studies) {
        studyEntries = new HashMap<>(studies.size());
        getImpl().setStudies(new ArrayList<>(studies.size()));
        for (StudyEntry study : studies) {
            getImpl().getStudies().add(study.getImpl());
            studyEntries.put(composeId(study.getStudyId()), study);
        }
    }

    @Deprecated
    public Map<String, StudyEntry> getSourceEntries() {
        return getStudiesMap();
    }

    public Map<String, StudyEntry> getStudiesMap() {
        if (getImpl().getStudies() != null) {
            if (studyEntries == null) {
                studyEntries = new HashMap<>();
                for (org.opencb.biodata.models.variant.avro.StudyEntry sourceEntry : getImpl().getStudies()) {
                    studyEntries.put(composeId(sourceEntry.getStudyId()), new StudyEntry(sourceEntry));
                }
            }
            return Collections.unmodifiableMap(studyEntries);
        }
        return null;
    }

    @Deprecated
    public StudyEntry getSourceEntry(String studyId) {
        return getStudy(studyId);
    }

    @Deprecated
    public StudyEntry getSourceEntry(String fileId, String studyId) {
        return getStudy(studyId);
    }

    public StudyEntry getStudy(String studyId) {
        if (getImpl().getStudies() != null) {
            return getStudiesMap().get(composeId(studyId));
        }
        return null;
    }

    public void addStudyEntry(StudyEntry studyEntry) {
        if (studyEntries == null) {
            studyEntries = new HashMap<>();
        }
        if (getImpl().getStudies() == null) {
            getImpl().setStudies(new ArrayList<>());
        }
        this.studyEntries.put(composeId(studyEntry.getStudyId()), studyEntry);
        getImpl().getStudies().add(studyEntry.getImpl());
    }

    public Iterable<String> getSampleNames(String studyId, String fileId) {
        StudyEntry file = getSourceEntry(studyId, fileId);
        if (file == null) {
            return null;
        }
        return file.getSamplesName();
    }

    public void transformToEnsemblFormat() {
        if (getType() == VariantType.INDEL || getType() == VariantType.SV || getLength() > 1) {
            if (getReference().charAt(0) == getAlternate().charAt(0)) {
                setReference(getReference().substring(1));
                setAlternate(getAlternate().substring(1));
                setStart(getStart() + 1);
                if (getReference().length() < getAlternate().length()) {
                    setEnd(getEnd() - 1);
                }

                if (getReference().equals("")) {
                    setReference("-");
                }
                if (getAlternate().equals("")) {
                    setAlternate("-");
                }

                resetLength();
            }
        }
    }


    @Override
    public String toString() {
        if (this.getReference() == null) {
            return getChromosome() + ":" + getStart() + ":" + (getAlternate().isEmpty() ? "-" : getAlternate());
        } else {
            return getChromosome() + ":" + getStart() + ":" + (getReference().isEmpty() ? "-" : getReference()) + ":" + (getAlternate().isEmpty() ? "-" : getAlternate());
        }
    }

    public String toJson() {
        return getImpl().toString();
    }

    public boolean sameGenomicVariant(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof Variant)) {
            return false;
        }

        Variant variant = (Variant) o;

        if (!Objects.equals(getStart(), variant.getStart())) {
            return false;
        }
        if (!Objects.equals(getEnd(), variant.getEnd())) {
            return false;
        }
        if (getChromosome() != null ? !getChromosome().equals(variant.getChromosome()) : variant.getChromosome() != null) {
            return false;
        }
        if (getReference() != null ? !getReference().equals(variant.getReference()) : variant.getReference() != null) {
            return false;
        }
        if (getAlternate() != null ? !getAlternate().equals(variant.getAlternate()) : variant.getAlternate() != null) {
            return false;
        }
        return getType() == variant.getType();

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Variant)) return false;

        Variant variant = (Variant) o;

        return !(getImpl() != null ? !getImpl().equals(variant.getImpl()) : variant.getImpl() != null);

    }

    @Override
    public int hashCode() {
        return getImpl() != null ? getImpl().hashCode() : 0;
    }


//    //    public int hashCode() {
//        int result = getChromosome() != null ? getChromosome().hashCode() : 0;
//        result = 31 * result + getStart();
//        result = 31 * result + getEnd();
//        result = 31 * result + (getReference() != null ? getReference().hashCode() : 0);
//        result = 31 * result + (getAlternate() != null ? getAlternate().hashCode() : 0);
//        result = 31 * result + (getType() != null ? getType().hashCode() : 0);
//        return result;
//    }

    private String composeId(String studyId) {
        return composeId(studyId, null);
    }

    @Deprecated
    private String composeId(String studyId, String fileId) {
        return studyId;
    }

    public static Variant parseVariant(String variantString) {
        return new Variant(variantString);
    }

    public static List<Variant> parseVariants(String variantsString) {
        List<Variant> variants = null;
        if(variantsString != null && !variantsString.isEmpty()) {
            String[] variantItems = variantsString.split(",");
            variants = new ArrayList<>(variantItems.length);
            for(String variantString: variantItems) {
                variants.add(new Variant(variantString));
            }
        }
        return variants;
    }

    public boolean overlapWith(Variant other, boolean inclusive) {
        return overlapWith(other.getChromosome(), other.getStart(), other.getEnd(), inclusive);
    }

    public boolean overlapWith(String chromosome, int start, int end, boolean inclusive) {
        if (!StringUtils.equals(this.getChromosome(), chromosome)) {
            return false; // Different Chromosome
        } else {
            int aStart = this.getStart();
            int aEnd = this.getEnd();

            if (aStart > aEnd) { // Insertion
                aStart = aEnd;
            }
            if (start > end){ // Insertion
                start = end;
            }

            if (inclusive) {
                return aStart <= end && aEnd >= start;
            } else {
                return aStart < end && aEnd > start;
            }
        }
    }

    public boolean onSameStartPosition (Variant other){
        return StringUtils.equals(this.getChromosome(), other.getChromosome()) 
                && this.getStart().equals(other.getStart());
    }

    /**
     * Check if Variant covers the same region (chromosome, start, end)
     * @param other Variant to check against
     * @return True if chromosome, start and end are the same
     */
    public boolean onSameRegion (Variant other){
        return onSameStartPosition(other) && this.getEnd().equals(other.getEnd());
    }

    /**
     * Return all VariantTypes subtypes given a VariantType.
     * {@link VariantType} represents a hierarchical structure where SNV includes all SNP, MNV includes MNP
     * and SV includes  INSERTION, DELETION, TRANSLOCATION, INVERSION and CNV
     *
     * @param variantType   Variant Type
     * @return  Set of subtypes
     */
    public static Set<VariantType> subTypes(VariantType variantType) {
        if(variantType.equals(VariantType.SNV)) {
            return Collections.singleton(VariantType.SNP);
        } else if (variantType.equals(VariantType.MNV)) {
            return Collections.singleton(VariantType.MNP);
        } else if (variantType.equals(VariantType.SV)) {
            return  SV_SUBTYPES;
        } else {
            return Collections.emptySet();
        }
    }

    @Override
    public int compareTo(Variant o) {
        if (this.equals(o)) {
            return 0;
        }
        return this.getImpl().compareTo(o.getImpl());
    }
}


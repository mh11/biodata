package org.opencb.biodata.tools.variant.converter;

import org.opencb.biodata.models.variant.StudyEntry;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.protobuf.VcfMeta;
import org.opencb.biodata.models.variant.protobuf.VcfSliceProtos;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Created on 05/11/15
 *
 * @author Jacobo Coll &lt;jacobo167@gmail.com&gt;
 */
public class VcfSliceToVariantListConverter implements Converter<VcfSliceProtos.VcfSlice, List<Variant>> {

    private final Map<String, Integer> samplesPosition;
    private final String fileId;
    private final String studyId;
    private final AtomicBoolean createMapCopy = new AtomicBoolean(false);

    public VcfSliceToVariantListConverter(VcfMeta meta) {
        this(meta.getVariantSource());
    }

    public VcfSliceToVariantListConverter(VariantSource source) {
        this(source.getSamplesPosition(), source.getFileId(), source.getStudyId());
    }

    public VcfSliceToVariantListConverter(Map<String, Integer> samplesPosition, String fileId, String studyId) {
//        this.samplesPosition = StudyEntry.sortSamplesPositionMap(samplesPosition);
        this.samplesPosition = new ConcurrentHashMap<>(samplesPosition);
        this.fileId = fileId;
        this.studyId = studyId;
    }

    public void setCreateMapCopy(boolean createMapCopy) {
        this.createMapCopy.set(createMapCopy);
    }

    @Override
    public List<Variant> convert(VcfSliceProtos.VcfSlice vcfSlice) {
        VcfRecordToVariantConverter recordConverter = new VcfRecordToVariantConverter(vcfSlice.getFields(),
                samplesPosition, fileId, studyId);
        recordConverter.setCreateMapCopy(isCreateMapCopy());
        List<Variant> variants = new ArrayList<>(vcfSlice.getRecordsCount());
        for (VcfSliceProtos.VcfRecord vcfRecord : vcfSlice.getRecordsList()) {
            variants.add(recordConverter.convert(vcfRecord, vcfSlice.getChromosome(), vcfSlice.getPosition()));
        }
        return variants;
    }

    private boolean isCreateMapCopy() {
        return this.createMapCopy.get();
    }

}

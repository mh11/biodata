package org.opencb.biodata.tools.variant.converter;

import com.google.protobuf.InvalidProtocolBufferException;
import org.apache.commons.lang3.StringUtils;
import org.junit.Before;
import org.junit.Test;
import org.opencb.biodata.models.variant.StudyEntry;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantNormalizer;
import org.opencb.biodata.models.variant.VariantVcfFactory;
import org.opencb.biodata.models.variant.avro.VariantType;
import org.opencb.biodata.models.variant.protobuf.VcfSliceProtos;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.opencb.biodata.models.variant.VariantTestUtils.generateVariantWithFormat;

/**
 * Created on 17/02/16
 *
 * @author Jacobo Coll &lt;jacobo167@gmail.com&gt;
 */
public class VcfSliceToVariantListConverterTest {


    private List<Variant> variants;

    @Before
    public void setUp() throws Exception {
        variants = Arrays.asList(
                generateVariantWithFormat("1:980:A:", "PASS", 102f,
                        toMap("K4", "V1", "K2", "V2", "END", "1000"), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1000:A:C", "PASS", 12f,
                        toMap("K3", "V1", "K4", "V2"), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1002:A:C", "PASS:LowGQX", 102f,
                        toMap("K5", "V1", "K2", "V2"), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1003:A:C", "PASS", 0f,
                        toMap("K3", "V1", "K2", "V2"), "GT:T", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1004:A:C", "LowGQX", null,
                        toMap("K2", "V1", "K3", "V2"), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1005:A:C", "PASS", 102f,
                        toMap("K3", "V1", "K2", "V2"), "GT:X", "S1", "0/1", "1"),
                generateVariantWithFormat("1:1006:A:", "PASS:LowGQX", 102f,
                        toMap("K1", "V1", "K5", "V2", "END", "1100"), "GT:T", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1106:T:C,TT", "PASS:LowGQX", 102f,
                        toMap("K2", "V1", "K3", "V2"), "GT:T", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1106:T:TT,C", "PASS:LowGQX", 102f,
                        toMap("K2", "V1", "K3", "V2"), "GT:T", "S1", "0/0", "1"),
                generateVariantWithFormat("1:1200:A:AATAG,AAGAAAGAG", "PASS", 102f,
                        toMap(), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:11111:A:", "PASS", 102f,
                        toMap(), "GT:X", "S1", "0/0", "1"),
                generateVariantWithFormat("1:11211::A,C", "PASS", 102f,
                        toMap(), "GT:X", "S1", "0/0", "1")
        );
        variants.get(5).getStudy("").getFile("").getAttributes().put(VariantVcfFactory.QUAL, ".");
        variants.get(6).setType(VariantType.NO_VARIATION);
        variants.get(6).setEnd(1100);

        variants.get(0).setType(VariantType.NO_VARIATION);
        variants.get(0).setEnd(1000);

        variants = new VariantNormalizer(true).apply(variants);
    }

    @Test
    public void buildFields() {

        VcfSliceProtos.Fields fields = VariantToVcfSliceConverter.buildDefaultFields(variants, null).build();

        assertEquals(Arrays.asList("PASS", "PASS:LowGQX", "LowGQX"), fields.getFiltersList());
        assertEquals(Arrays.asList("GT:X", "GT:T"), fields.getFormatsList());
        assertEquals(Arrays.asList("K2", "K3", "K4", "K5", "K1"), fields.getInfoKeysList());
        assertEquals(Arrays.asList(0, 1), fields.getDefaultInfoKeysList());
        assertEquals(Arrays.asList("0/0", "0/1"), fields.getGtsList());

    }

    @Test
    public void testParallelConversion() throws InvalidProtocolBufferException {

        LinkedHashMap<String, Integer> samplesPosition = variants.get(0).getStudies().get(0).getSamplesPosition();
        VcfSliceToVariantListConverter vcfSliceToVariantListConverter = new VcfSliceToVariantListConverter(samplesPosition, "", "");

        VariantToVcfSliceConverter converter = new VariantToVcfSliceConverter();
        final byte[] sliceBytes = converter.convert(variants, 1000).toByteArray();

        IntStream.range(0, 100).parallel().forEach(x -> {
            List<Variant> collect = IntStream.range(0, 100).parallel().boxed().flatMap(i -> {
                try {
                    VcfSliceProtos.VcfSlice slice = VcfSliceProtos.VcfSlice.parseFrom(sliceBytes);
                    List<Variant> convert = vcfSliceToVariantListConverter.convert(slice);

                    assertEquals(VariantType.NO_VARIATION, convert.get(6).getType());
                    assertEquals(1000, convert.get(0).getEnd().intValue());
                    assertEquals(1000, convert.get(1).getEnd().intValue());
                    assertEquals(1100, convert.get(6).getEnd().intValue());
                    assertEquals("0", convert.get(3).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));
                    assertEquals(null, convert.get(4).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));
                    assertEquals(null, convert.get(5).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));
                    return convert.stream();
                } catch (InvalidProtocolBufferException e) {
                    throw new IllegalStateException(e);
                }
            }).collect(Collectors.toList());

            collect.parallelStream().forEach(v -> {
                Map<String, String> gt = sampleToSampleData(v, "GT");
                assertEquals(gt.size(), 1);
                assertNotNull(gt.get("S1"));
                StudyEntry studyEntry = getStudy(v);
                String sampleData = studyEntry.getSampleData("S1", "GT");
                assertNotNull("Sampel GT are null", sampleData);
            });

        });
    }

    private Map<String, String> sampleToSampleData(Variant var, String key){
        StudyEntry se = getStudy(var);
        return se.getSamplesName().stream()
                .filter(e -> StringUtils.isNotBlank(se.getSampleData(e, key))) // check for NULL or empty string
                .collect(Collectors.toMap(e -> e, e -> se.getSampleData(e, key)));
    }

    static StudyEntry getStudy(Variant variant) {
        return variant.getStudies().get(0);
    }

    @Test
    public void testConvertVariants() throws InvalidProtocolBufferException {
//        VcfSliceToVariantListConverter converter = new VcfSliceToVariantListConverter();
        VariantToVcfSliceConverter converter = new VariantToVcfSliceConverter();
        VcfSliceProtos.VcfSlice slice = converter.convert(variants, 1000);

        slice = VcfSliceProtos.VcfSlice.parseFrom(slice.toByteArray());

        LinkedHashMap<String, Integer> samplesPosition = variants.get(0).getStudies().get(0).getSamplesPosition();
        VcfSliceToVariantListConverter vcfSliceToVariantListConverter = new VcfSliceToVariantListConverter(samplesPosition, "", "");
        List<Variant> convert = vcfSliceToVariantListConverter.convert(slice);


        assertEquals(VariantType.NO_VARIATION, convert.get(6).getType());
        assertEquals(1000, convert.get(0).getEnd().intValue());
        assertEquals(1000, convert.get(1).getEnd().intValue());
        assertEquals(1100, convert.get(6).getEnd().intValue());
        assertEquals("0", convert.get(3).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));
        assertEquals(null, convert.get(4).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));
        assertEquals(null, convert.get(5).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL));

        for (int i = 0; i < convert.size(); i++) {
            // Set qual to NULL if required.
            if (".".equals(variants.get(i).getStudy("").getFile("").getAttributes().get(VariantVcfFactory.QUAL))) {
                variants.get(i).getStudy("").getFile("").getAttributes().put(VariantVcfFactory.QUAL, null);
            }

            System.out.println("Expected  : " + variants.get(i).toJson());
            System.out.println("Converted : " + convert.get(i).toJson());
            System.out.println("Proto     : " + slice.getRecords(i).toString());
            assertEquals(variants.get(i), convert.get(i));
            System.out.println("------------------------------");
        }
    }

    public static Map<String, String> toMap(String... strings) {
        assertEquals(0, strings.length % 2);
        Map<String, String> map = new HashMap<>(strings.length / 2);
        for (int i = 0; i < strings.length; i = i + 2) {
            map.put(strings[i], strings[i + 1]);
        }
        return map;
    }

}
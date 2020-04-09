package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.DuplicateSetWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import picard.cmdline.programgroups.OtherProgramGroup;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = OtherProgramGroup.class
)
public class CorrectErrorsInDuplicateReads extends DuplicateSetWalker {
    @Override
    public void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext) {
    }
}

package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;

public interface ErrorCorrector {
    public Collection<GATKRead> correctErrorsAndUpdateQualities(final DuplicateSet duplicateSet, final ReferenceContext referenceContext);

    public boolean filterDuplicateSet(final DuplicateSet duplicateSet);
}

package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AdaptiveChainPruner<V extends BaseVertex, E extends BaseEdge> extends ChainPruner<V,E> {
    private final double initialErrorProbability;
    private final double logOddsThreshold;
    private final int maxUnprunedVariants;

    public AdaptiveChainPruner(final double initialErrorProbability, final double logOddsThreshold, final int maxUnprunedVariants) {
        ParamUtils.isPositive(initialErrorProbability, "Must have positive error probability");
        this.initialErrorProbability = initialErrorProbability;
        this.logOddsThreshold = logOddsThreshold;
        this.maxUnprunedVariants = maxUnprunedVariants;
    }

    @Override
    protected Collection<Path<V,E>> chainsToRemove(final List<Path<V, E>> chains) {
        if (chains.isEmpty()) {
            return Collections.emptyList();
        }

        final BaseGraph<V,E> graph = chains.get(0).getGraph();



        Collection<Path<V,E>> probableErrorChains = likelyErrorChains(chains, graph, initialErrorProbability);
        final int errorCount = probableErrorChains.stream().mapToInt(c -> c.getLastEdge().getMultiplicity()).sum();
        final int totalBases = chains.stream().mapToInt(c -> c.getEdges().stream().mapToInt(E::getMultiplicity).sum()).sum();
        final double errorRate = (double) errorCount / totalBases;

        return likelyErrorChains(chains, graph, errorRate).stream().filter(c -> !c.getEdges().stream().anyMatch(BaseEdge::isRef)).collect(Collectors.toList());
    }

    private Collection<Path<V,E>> likelyErrorChains(final List<Path<V, E>> chains, final BaseGraph<V,E> graph, final double errorRate) {

        // start from set of really obvious good chains, defined as those with an edge with multiplicity above some fraction of the
        // maximum multiplicity edge.  Then grow the set of good chains by looking at the log odds of junctions between
        // good chains and unknown chains.
        final int maxEdgeWeight = graph.edgeSet().stream().mapToInt(BaseEdge::getMultiplicity).max().orElse(0);
        final int thresholdEdgeWeight = maxEdgeWeight / 10;
        final Set<Path<V, E>> definiteGoodChains = chains.stream()
                .filter(chain -> chain.getEdges().stream().mapToInt(BaseEdge::getMultiplicity).max().orElse(0) > thresholdEdgeWeight)
                .collect(Collectors.toSet());

        final Set<V> terminalVerticesOfGoodChains = definiteGoodChains.stream()
                .flatMap(c -> Stream.of(c.getFirstVertex(), c.getLastVertex())).collect(Collectors.toSet());

        // initialize result to all chains that aren't definitely good and iteratively reduce this set as chains in
        // question are found to have good log odds at junctions with good chains
        final Set<Path<V,E>> errorChains = chains.stream().filter(c -> !definiteGoodChains.contains(c)).collect(Collectors.toSet());

        // we only need to calculate log odds for questionable chains
        final Map<Path<V,E>, Pair<Double, Double>> chainLogOdds = errorChains.stream()
                .collect(Collectors.toMap(c -> c, c-> chainLogOdds(c, graph, errorRate)));

        while (true) {
            final Set<Path<V, E>> moreGoodChains = new HashSet<>();

            for (final Path<V, E> chain : errorChains) {
                if (terminalVerticesOfGoodChains.contains(chain.getFirstVertex()) && chainLogOdds.get(chain).getLeft() > logOddsThreshold) {
                    moreGoodChains.add(chain);
                    terminalVerticesOfGoodChains.add(chain.getLastVertex());
                } else if (terminalVerticesOfGoodChains.contains(chain.getLastVertex()) && chainLogOdds.get(chain).getRight() > logOddsThreshold) {
                    moreGoodChains.add(chain);
                    terminalVerticesOfGoodChains.add(chain.getFirstVertex());
                }
            }

            if (moreGoodChains.isEmpty()) {
                break;
            } else {
                errorChains.removeAll(moreGoodChains);
            }
        }

        // add non-error chains to error chains if maximum number of variants has been exceeded
        chains.stream()
                .filter(c -> !errorChains.contains(c))
                .filter(c -> !definiteGoodChains.contains(c))
                .filter(c -> isChainPossibleVariant(c, graph))
                .sorted(Comparator.comparingDouble(c -> Math.min(chainLogOdds.get(c).getLeft(), chainLogOdds.get(c).getLeft())).reversed())
                .skip(maxUnprunedVariants)
                .forEach(errorChains::add);

        return errorChains;

    }

    // left and right chain log odds
    private Pair<Double, Double> chainLogOdds(final Path<V,E> chain, final BaseGraph<V,E> graph, final double errorRate) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        final double leftLogOdds = graph.isSource(chain.getFirstVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
        final double rightLogOdds = graph.isSink(chain.getLastVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

        return ImmutablePair.of(leftLogOdds, rightLogOdds);
    }
    
    private boolean isChainPossibleVariant(final Path<V,E> chain, final BaseGraph<V,E> graph) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        return leftMultiplicity <= leftTotalMultiplicity / 2 || rightMultiplicity <= rightTotalMultiplicity / 2;
    }
}

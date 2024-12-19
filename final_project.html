<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en"><head>

<meta charset="utf-8">
<meta name="generator" content="quarto-1.5.57">

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">


<title>final_project</title>
<style>
code{white-space: pre-wrap;}
span.smallcaps{font-variant: small-caps;}
div.columns{display: flex; gap: min(4vw, 1.5em);}
div.column{flex: auto; overflow-x: auto;}
div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
ul.task-list{list-style: none;}
ul.task-list li input[type="checkbox"] {
  width: 0.8em;
  margin: 0 0.8em 0.2em -1em; /* quarto-specific, see https://github.com/quarto-dev/quarto-cli/issues/4556 */ 
  vertical-align: middle;
}
/* CSS for syntax highlighting */
pre > code.sourceCode { white-space: pre; position: relative; }
pre > code.sourceCode > span { line-height: 1.25; }
pre > code.sourceCode > span:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode > span { color: inherit; text-decoration: inherit; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
pre > code.sourceCode { white-space: pre-wrap; }
pre > code.sourceCode > span { display: inline-block; text-indent: -5em; padding-left: 5em; }
}
pre.numberSource code
  { counter-reset: source-line 0; }
pre.numberSource code > span
  { position: relative; left: -4em; counter-increment: source-line; }
pre.numberSource code > span > a:first-child::before
  { content: counter(source-line);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
  }
pre.numberSource { margin-left: 3em;  padding-left: 4px; }
div.sourceCode
  {   }
@media screen {
pre > code.sourceCode > span > a:first-child::before { text-decoration: underline; }
}
</style>


<script src="final_project_files/libs/clipboard/clipboard.min.js"></script>
<script src="final_project_files/libs/quarto-html/quarto.js"></script>
<script src="final_project_files/libs/quarto-html/popper.min.js"></script>
<script src="final_project_files/libs/quarto-html/tippy.umd.min.js"></script>
<script src="final_project_files/libs/quarto-html/anchor.min.js"></script>
<link href="final_project_files/libs/quarto-html/tippy.css" rel="stylesheet">
<link href="final_project_files/libs/quarto-html/quarto-syntax-highlighting.css" rel="stylesheet" id="quarto-text-highlighting-styles">
<script src="final_project_files/libs/bootstrap/bootstrap.min.js"></script>
<link href="final_project_files/libs/bootstrap/bootstrap-icons.css" rel="stylesheet">
<link href="final_project_files/libs/bootstrap/bootstrap.min.css" rel="stylesheet" id="quarto-bootstrap" data-mode="light">


</head>

<body class="fullcontent">

<div id="quarto-content" class="page-columns page-rows-contents page-layout-article">

<main class="content" id="quarto-document-content">




<section id="robustness-of-phylogenetic-tree-inference---apoditrysia" class="level1">
<h1>Robustness of Phylogenetic Tree Inference - Apoditrysia</h1>
</section>
<section id="introduction-and-goals" class="level1">
<h1>Introduction and Goals</h1>
<p>This project aims to examine the sensitivity of the phylogenetic trees inferred from Apoditrysia (butterflies and larger moths) sequences. Two datasets were created. The first dataset consisted of 18S sequences from Apoditrysia, which are highly conserved and therefore useful for phylogenetic analyses. The second dataset contained full mitochondrial genomes. Both datasets were extracted from BLAST using a bait sequence and aligned using MAFFT.</p>
<p>Dataset 1:<br>
To begin, 166 sequences were found in Apoditrysia. The outgroup taxa were chosen such that 6/8 occured in Ditrysia clade but not the Apoditrysia clade. Then the remaining two outgroup taxa were from Incurvarioidea (yucca and fairy moths), which is in Heteroneura along with Ditrysia. The goal of such outgroup sampling was to find outgroups that have varying degrees of similarity to the target clade.</p>
<p>In summary<br>
Outgroup: 6/8 in Ditrysia.<br>
Gelechioidea: Phthorimaea_operculella<br>
Obtectomera: Loryma_sp.<br>
Tineoidea: Ogmograptis_sp., Thyridopteryx_ephemeraeformis<br>
Yponomeutoidea: Tritymba_sp., Plutella_xylostella<br>
</p>
<p>2/8 in sister taxa to Incurvarioidea<br>
Prodoxus_quinquepunctellus, Tegeticula_yuccasella</p>
<p>Dataset 2:<br>
96 species from Apoditrysia<br>
5 outgroup species: Polyploca_ridens, Platyedra_subcinerea, Parnassius_tianschanicus, Plutella_armoraciae, Paraclemensia_caerulea from groups with various distances from Apoditrysia.<br>
</p>
<p>The sequences were then combined, and regex was used to adjust the labels of the FASTA files.</p>
<p>The particular regex used was Find: &gt;\S* (\S<em>)\s(\S</em>).* Replace: &gt;$1_$2</p>
<p>Then the sequences were aligned using MAFFT and analyses were performed using the McCleary Cluster.</p>
</section>
<section id="methods" class="level1">
<h1>Methods</h1>
<p>IQTree and RevBayes were used to generate trees. The R package TreeDist was used to quantify and visualize the difference between trees.</p>
<p>All input fasta files, .sh files, .rev files, and output files are included in this directory.</p>
</section>
<section id="results" class="level1">
<h1>Results</h1>
<p>IQTree was first used to generate the most likely tree using bootstrapping (1000 iterations) on all extracted 18S sequences.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb1"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a>first_pass <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"all_18s_aligned.fa.treefile"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb2"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(first_pass)</span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(first_pass<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-2-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>It may not be very obvious, but the bootstrap support of this tree is not great. Below is a histogram of the supports.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb3"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(first_pass<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-3-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<p>Therefore, the dataset was pruned by removing sequences that were much longer (i.e.&nbsp;3-4k nucleotides longer than the majority). This was about 10 sequences and helped the alignment to have less gaps.</p>
<p>(Previously, most alignmemts had ~60% gaps on average).</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb4"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb4-1"><a href="#cb4-1" aria-hidden="true" tabindex="-1"></a>pruned <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"18s_pruned_aligned.fasta.treefile"</span>)</span>
<span id="cb4-2"><a href="#cb4-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(pruned)</span>
<span id="cb4-3"><a href="#cb4-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(pruned<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-4-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>And a histogram of the bootstrap supports for this tree:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb5"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(pruned<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-5-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<p>As one can see, while the support it is still not ideal, the histogram appears to be much better.</p>
<p>Now, to do another sanity check, one can observe the location of the outgroup taxa (in red).</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb6"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a>colors <span class="ot">=</span> <span class="fu">ifelse</span>(<span class="fu">grepl</span>(<span class="st">"Phthorimaea_operculella|Loryma_sp.|Ogmograptis_sp.|Thyridopteryx_ephemeraeformis|Tritymba_sp.|Plutella_xylostella|Prodoxus_quinquepunctellus|Tegeticula_yuccasella"</span>,pruned<span class="sc">$</span>tip.label),<span class="st">"red"</span>,<span class="st">"black"</span>)</span>
<span id="cb6-2"><a href="#cb6-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-3"><a href="#cb6-3" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(pruned,<span class="at">tip.color=</span>colors)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-6-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Most of the outgroup taxa are near each other as desired. However, they are in the middle of the phylogeny for some reason (perhaps this is because MAFFT adjusted the order). However, while this “looks” weird, it is not actually significant to the topology of the tree, so I will reorder for future analysis.</p>
<p>More concerning is that <em>Thyridopteryx</em> is clearly far removed from the remainder of its outgroup brethren. Further inverstigation reveals that it is not in Tineoidea as BLAST promised; therefore, a more closely related species to Apoditrysia might be chosen instead.</p>
<p>Finally, an analysis of the sequences reveals that <em>Epichoristodes_acerbella</em> has a very sparce alignment. Thus, I will remove it as it is uninformative in its current state. This can be seen in the produced tree in which Epichoristodes stands out from the rest of the taxa.</p>
<p>Below is the updated tree with <em>Epichoristodes_acerbella</em> removed. The bootstrap support is better but still not ideal, and the outgroups are generally together.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb7"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb7-1"><a href="#cb7-1" aria-hidden="true" tabindex="-1"></a>pruned2 <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"v2_18s_pruned_aligned.fasta.treefile"</span>)</span>
<span id="cb7-2"><a href="#cb7-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(pruned2)</span>
<span id="cb7-3"><a href="#cb7-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(pruned2<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-7-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb8"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb8-1"><a href="#cb8-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(pruned2<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-8-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb9"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb9-1"><a href="#cb9-1" aria-hidden="true" tabindex="-1"></a>colors <span class="ot">=</span> <span class="fu">ifelse</span>(<span class="fu">grepl</span>(<span class="st">"Phthorimaea_operculella|Loryma_sp.|Ogmograptis_sp.|Thyridopteryx_ephemeraeformis|Tritymba_sp.|Plutella_xylostella|Prodoxus_quinquepunctellus|Tegeticula_yuccasella"</span>,pruned2<span class="sc">$</span>tip.label),<span class="st">"red"</span>,<span class="st">"black"</span>)</span>
<span id="cb9-2"><a href="#cb9-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb9-3"><a href="#cb9-3" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(pruned2,<span class="at">tip.color=</span>colors)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-9-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Additionally, IQTree has an outgroup function that should (theoretically) group all the outgroups together. However, the tree generated using this method not only has lower bootstrap support values but also does not group the outgroup “clade” together. Therefore, the previous model is condsidered preferable.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb10"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb10-1"><a href="#cb10-1" aria-hidden="true" tabindex="-1"></a>t <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"outgroup_auto.fasta.treefile"</span>)</span>
<span id="cb10-2"><a href="#cb10-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(t)</span>
<span id="cb10-3"><a href="#cb10-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(t<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-10-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb11"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb11-1"><a href="#cb11-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(t<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-11-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb12"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb12-1"><a href="#cb12-1" aria-hidden="true" tabindex="-1"></a>colors <span class="ot">=</span> <span class="fu">ifelse</span>(<span class="fu">grepl</span>(<span class="st">"Phthorimaea_operculella|Loryma_sp.|Ogmograptis_sp.|Thyridopteryx_ephemeraeformis|Tritymba_sp.|Plutella_xylostella|Prodoxus_quinquepunctellus|Tegeticula_yuccasella"</span>,t<span class="sc">$</span>tip.label),<span class="st">"red"</span>,<span class="st">"black"</span>)</span>
<span id="cb12-2"><a href="#cb12-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb12-3"><a href="#cb12-3" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(t,<span class="at">tip.color=</span>colors)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-12-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Next, the preferred ML tree was compared to a Bayesian tree created using RevBayes (HKY).</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb13"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb13-1"><a href="#cb13-1" aria-hidden="true" tabindex="-1"></a>t <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"bayes_HKY.nex"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="bayes_HKY.tre.jpg" class="class figure-img" style="width:100.0%;height:125.0%"></p>
<figcaption>HKY 18S</figcaption>
</figure>
</div>
<p>Note that the posterior probability of each branch is very low – most of them are below 0.2.</p>
<p>Therefore, it remained clear that more data was needed in order to create a more robust tree. Espeland et al.&nbsp;(2018) created a phylogeny of butterflies (closely related to Apoditrysia) and provided a list of probes. Python was used to randomly sample 15% of the probes. The sampled probes were then used in nBLAST. It was observed that the CO1 gene (Cytochrome c oxidase I) was well-sampled within <em>Apoditrysia</em>; moreover, a number of species had complete mitochondria gene sequences that contained CO1. Therefore, a complete mitochondria sequence was used as a probe in nBLAST.</p>
<p>The resulting IQtree model had much better results: the bootstrap supports were higher. However, the outgroups were not next to each other; this could be a result of choosing outgroups too close to the ingroup (in fact, NCBI does not consider those in Obtectomera to be within Apoditrysia but Wikipedia does.)</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb14"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb14-1"><a href="#cb14-1" aria-hidden="true" tabindex="-1"></a>mito <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"mito1_al.fasta.treefile"</span>)</span>
<span id="cb14-2"><a href="#cb14-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(mito)</span>
<span id="cb14-3"><a href="#cb14-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(mito<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-14-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb15"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb15-1"><a href="#cb15-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(mito<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-15-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb16"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb16-1"><a href="#cb16-1" aria-hidden="true" tabindex="-1"></a>colors <span class="ot">=</span> <span class="fu">ifelse</span>(<span class="fu">grepl</span>(<span class="st">"Polyploca_ridens|Platyedra_subcinerea|Parnassius_tianschanicus|Plutella_armoraciae|Paraclemensia_caerulea"</span>,mito<span class="sc">$</span>tip.label),<span class="st">"red"</span>,<span class="st">"black"</span>)</span>
<span id="cb16-2"><a href="#cb16-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb16-3"><a href="#cb16-3" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(mito,<span class="at">tip.color=</span>colors)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-16-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Thus another IQTree run was completed to see if the outgroups had any impact on the tree at all. Below is the no-outgroup tree which apppears very similar to the above standard auto tree for the mitochondria data.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb17"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb17-1"><a href="#cb17-1" aria-hidden="true" tabindex="-1"></a>x_out <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"mito1_nooutgroup.fasta.treefile"</span>)</span>
<span id="cb17-2"><a href="#cb17-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(x_out)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-17-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>The TreeDist package was then used to quantify the differences between these 2 trees. In particular, the TreeDistance function shows the clustering information variation between the mitochondial trees with and without outgroups normalized by the total information content of all splits.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb18"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb18-1"><a href="#cb18-1" aria-hidden="true" tabindex="-1"></a>out_diff <span class="ot">=</span> <span class="fu">TreeDistance</span>(mito, x_out)</span>
<span id="cb18-2"><a href="#cb18-2" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(out_diff)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>[1] 0.03208044</code></pre>
</div>
</div>
<p>which indicates that the trees are very similar.</p>
<p>Visually, the difference between the trees was as follows:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb20"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb20-1"><a href="#cb20-1" aria-hidden="true" tabindex="-1"></a><span class="fu">VisualizeMatching</span>(MatchingSplitDistance, mito, x_out, <span class="at">edge.cex=</span><span class="cn">FALSE</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stderr">
<pre><code>Warning in edge.width[se] &lt;- 1 + (10 * ns): number of items to replace is not a
multiple of replacement length</code></pre>
</div>
<div class="cell-output cell-output-stderr">
<pre><code>Warning in edge.color[se] &lt;- edgeColPalette[1 + ceiling(255 * ns)]: number of
items to replace is not a multiple of replacement length</code></pre>
</div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-19-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>In the above example, the yellow colored splits indicate a difference between both trees. This only occurs due to the removal of the outgroups, implying that the outgroups have little to no impact on the tree itself.</p>
<p>So how sensitive is the is the base ML tree to ablation? 25 species were removed at random from the tree, and two ways of re-aligning the data were completed. The first was a “clipping” of previous alignments; that is, the alignment of the 75 remaining species did not change, and less sequences were used to construct the tree. The second is an actual re-alignment of all kept species.</p>
<p>Here is the clipped version:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb23"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb23-1"><a href="#cb23-1" aria-hidden="true" tabindex="-1"></a>clipped <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"mito_abl_clipped.fasta.treefile"</span>)</span>
<span id="cb23-2"><a href="#cb23-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(clipped)</span>
<span id="cb23-3"><a href="#cb23-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(clipped<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-20-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb24"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb24-1"><a href="#cb24-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(clipped<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-21-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<p>And here is the realigned version:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb25"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb25-1"><a href="#cb25-1" aria-hidden="true" tabindex="-1"></a>realigned <span class="ot">=</span> <span class="fu">read.tree</span>(<span class="st">"mito_abl.fasta.treefile"</span>)</span>
<span id="cb25-2"><a href="#cb25-2" aria-hidden="true" tabindex="-1"></a><span class="fu">plot</span>(realigned)</span>
<span id="cb25-3"><a href="#cb25-3" aria-hidden="true" tabindex="-1"></a><span class="fu">nodelabels</span>(realigned<span class="sc">$</span>node.label)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-22-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb26"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb26-1"><a href="#cb26-1" aria-hidden="true" tabindex="-1"></a><span class="fu">hist</span>(<span class="fu">as.integer</span>(realigned<span class="sc">$</span>node.label), <span class="at">xlab=</span><span class="st">"Bootstrap support"</span>, <span class="at">main=</span><span class="st">""</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-23-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<p>Now we can compare the trees.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb27"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb27-1"><a href="#cb27-1" aria-hidden="true" tabindex="-1"></a>tocomp_trees <span class="ot">=</span> <span class="fu">structure</span>(<span class="fu">list</span>(<span class="at">base =</span> mito, <span class="at">clipped =</span> clipped , <span class="at">realigned =</span> realigned), <span class="at">class =</span> <span class="st">'multiPhylo'</span>)</span>
<span id="cb27-2"><a href="#cb27-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb27-3"><a href="#cb27-3" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(<span class="fu">TreeDistance</span>(tocomp_trees, tocomp_trees))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>                base    clipped realigned
base      0.00000000 0.08843369 0.2960598
clipped   0.06401343 0.00000000 0.2378496
realigned 0.27654539 0.23715756 0.0000000</code></pre>
</div>
</div>
<p>Clearly, the realigned version has more differences to the base and clipped trees than the base tree has to the clipped tree. This indicates that IQTree is very sensitive to the alignment of the data.</p>
<p>We can also visually compare the clipped and realigned trees using the same distance metric as the previous example (that is, recording the split differences of the two trees)</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb29"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb29-1"><a href="#cb29-1" aria-hidden="true" tabindex="-1"></a><span class="fu">VisualizeMatching</span>(MatchingSplitDistance, clipped, realigned, <span class="at">edge.cex=</span><span class="cn">FALSE</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-25-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>While there are similarities between the two trees, there are also some differences especially pertaining to the upper part of the tree.</p>
<p>It is also interesting to compare the IQTrees to those created using RevBayes. First, two trees were created using the ‘base’ dataset according to the HKY and GTR models.</p>
<p>HKY:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb30"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb30-1"><a href="#cb30-1" aria-hidden="true" tabindex="-1"></a>rb_hky <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"rb_mito1_al_HKY.tre"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="rb_mito1_al_HKY.tre.jpg" class="class figure-img" style="width:100.0%;height:110.0%"></p>
<figcaption>HKY Base Model</figcaption>
</figure>
</div>
<p>Note that similar to the IQTree trees, this HKY model shows vast improvements to the posterior probabilities compared to the 18S data.</p>
<p>GTR:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb31"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb31-1"><a href="#cb31-1" aria-hidden="true" tabindex="-1"></a>rb_gtr <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"rb_mito1_al_GTR.tre"</span>)</span>
<span id="cb31-2"><a href="#cb31-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb31-3"><a href="#cb31-3" aria-hidden="true" tabindex="-1"></a><span class="co"># removing single quotes</span></span>
<span id="cb31-4"><a href="#cb31-4" aria-hidden="true" tabindex="-1"></a><span class="cf">for</span> (i <span class="cf">in</span> <span class="dv">1</span><span class="sc">:</span><span class="fu">length</span>(rb_gtr<span class="sc">$</span>tip.label)) {</span>
<span id="cb31-5"><a href="#cb31-5" aria-hidden="true" tabindex="-1"></a>  s <span class="ot">&lt;-</span> rb_gtr<span class="sc">$</span>tip.label[i]</span>
<span id="cb31-6"><a href="#cb31-6" aria-hidden="true" tabindex="-1"></a>  <span class="cf">if</span> (<span class="fu">substring</span>(s, <span class="dv">1</span>, <span class="dv">1</span>) <span class="sc">==</span> <span class="st">"'"</span>) {</span>
<span id="cb31-7"><a href="#cb31-7" aria-hidden="true" tabindex="-1"></a>    rb_gtr<span class="sc">$</span>tip.label[i] <span class="ot">&lt;-</span> <span class="fu">substring</span>(s, <span class="dv">2</span>, <span class="fu">nchar</span>(s) <span class="sc">-</span> <span class="dv">1</span>)</span>
<span id="cb31-8"><a href="#cb31-8" aria-hidden="true" tabindex="-1"></a>  }</span>
<span id="cb31-9"><a href="#cb31-9" aria-hidden="true" tabindex="-1"></a>}</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="rb_mito1_al_GTR.tre.jpg" class="class figure-img" style="width:100.0%;height:120.0%"></p>
<figcaption>GTR Base BI Model</figcaption>
</figure>
</div>
<p>Comparing these Bayesian-generated base trees to the ML-generated IQTrees from before:</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb32"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb32-1"><a href="#cb32-1" aria-hidden="true" tabindex="-1"></a>tocomp_trees <span class="ot">=</span> <span class="fu">structure</span>(<span class="fu">list</span>(<span class="at">ML_Base =</span> mito, <span class="at">RB_HKY =</span> rb_hky , <span class="at">RB_GTR =</span> rb_gtr), <span class="at">class =</span> <span class="st">'multiPhylo'</span>)</span>
<span id="cb32-2"><a href="#cb32-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb32-3"><a href="#cb32-3" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(<span class="fu">TreeDistance</span>(tocomp_trees, tocomp_trees))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>          ML_Base     RB_HKY     RB_GTR
ML_Base 0.0000000 0.18513679 0.18767012
RB_HKY  0.1851368 0.00000000 0.09500635
RB_GTR  0.1876701 0.09500635 0.00000000</code></pre>
</div>
</div>
<p>Note that the metric we are using here (Information-based generalized Robinson–Foulds distances) is not symmetric; the distance from ML_Base to RB_HKY is not equal to the distance from RB_HKY to ML_Base. Clearly, RF distance is symmetric, so the cause of this asymmetry is probably a difference in the information content of some splits in the trees.</p>
<p>However, the general pattern remains clear. The RevBayes trees are closer to each other than to the IQTree. This suggests that this particular dataset could be slightly sensitive to the phylogenetic tree-generating method.</p>
<p>Next, the two RevBayes base trees were visually compared to analyze where their differences occur. It is clear that most species are in similar places in both trees (for example, <em>Pidorus</em> moved within the <em>Pidorus</em>, <em>Histia</em>, <em>Eterusia</em>, <em>Erasmia</em>, <em>Amesia</em> clade).</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb34"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb34-1"><a href="#cb34-1" aria-hidden="true" tabindex="-1"></a><span class="fu">VisualizeMatching</span>(MatchingSplitDistance, rb_hky, rb_gtr, <span class="at">edge.cex=</span><span class="cn">FALSE</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-29-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Note that <em>Grapholita</em> is not monophyletic in the HKY tree but is monophyletic in the GTR tree. Given that the latter is probably more likely, further ablation and no-outgroup trees were built using RevBayes GTR.</p>
<p>Below is the no outgroup Bayesian Infererance GTR tree.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb35"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb35-1"><a href="#cb35-1" aria-hidden="true" tabindex="-1"></a>rb_xog <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"rb_mito1_nooutgroup_GTR.tre"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="rb_mito1_nooutgroup_GTR.tre.jpg" class="class figure-img" style="width:100.0%;height:125.0%"></p>
<figcaption>GTR No Outgroup BI Model</figcaption>
</figure>
</div>
<p>Notably, removing the outgroup lowers the posterior probabilities of many splits. However, the resulting trees are very similar (a distance of ~0.08 as seen below).</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb36"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb36-1"><a href="#cb36-1" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(<span class="fu">TreeDistance</span>(rb_gtr, rb_xog))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>[1] 0.08107694</code></pre>
</div>
</div>
<p>Then the two ablation trees (clipped and realignment) were built.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb38"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb38-1"><a href="#cb38-1" aria-hidden="true" tabindex="-1"></a>rb_clip <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"rb_mito_abl_clipped_GTR.tre"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="rb_mito_abl_clipped_GTR.tre.jpg" class="class figure-img" style="width:100.0%;height:125.0%"></p>
<figcaption>GTR Clipped Ablation BI Model</figcaption>
</figure>
</div>
<div class="cell">
<div class="sourceCode cell-code" id="cb39"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb39-1"><a href="#cb39-1" aria-hidden="true" tabindex="-1"></a>rb_realign <span class="ot">=</span> <span class="fu">read.nexus</span>(<span class="st">"rb_mito_abl_GTR.tre"</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</div>
<div id="id" class="quarto-figure quarto-figure-center anchored">
<figure class="figure">
<p><img src="rb_mito_abl_GTR.tre.jpg" class="class figure-img" style="width:100.0%;height:125.0%"></p>
<figcaption>GTR Realigned Ablation BI Model</figcaption>
</figure>
</div>
<p>Clearly, both have slightly reduced posterior probablilities compared to the base BI model.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb40"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb40-1"><a href="#cb40-1" aria-hidden="true" tabindex="-1"></a>tocomp_trees <span class="ot">=</span> <span class="fu">structure</span>(<span class="fu">list</span>(<span class="at">BI_Base =</span> rb_gtr, <span class="at">Clipped =</span> rb_clip , <span class="at">Realigned =</span> rb_realign), <span class="at">class =</span> <span class="st">'multiPhylo'</span>)</span>
<span id="cb40-2"><a href="#cb40-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb40-3"><a href="#cb40-3" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(<span class="fu">TreeDistance</span>(tocomp_trees, tocomp_trees))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output cell-output-stdout">
<pre><code>            BI_Base    Clipped  Realigned
BI_Base   0.0000000 0.08856116 0.13683887
Clipped   0.1014466 0.00000000 0.07444176
Realigned 0.1377383 0.06214738 0.00000000</code></pre>
</div>
</div>
<p>However, the clipped and realigned models are very close to one another (and both are closer to BI base GTR than either the ML clipped or ML realigned were to the base ML model). This suggests that perhaps Bayesian inference is a better fit for this data set as the model is more robust.</p>
<p>Finally, the clipped and realigned models are visualized to contextualize their differences.</p>
<div class="cell">
<div class="sourceCode cell-code" id="cb42"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb42-1"><a href="#cb42-1" aria-hidden="true" tabindex="-1"></a><span class="fu">VisualizeMatching</span>(MatchingSplitDistance, rb_clip, rb_realign, <span class="at">edge.cex=</span><span class="cn">FALSE</span>)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="final_project_files/figure-html/unnamed-chunk-35-1.png" class="img-fluid figure-img" width="864"></p>
</figure>
</div>
</div>
</div>
<p>Again, one observes that the trees are very similar with small clade rearrangements.</p>
</section>
<section id="discussion" class="level1">
<h1>Discussion</h1>
<p>In summary, the quality (and quantity) of the data used in any phylogenetic inference model is vital to the overall robustness of the model. More than one gene is therefore preferable to build a strong model. This was shown via the difference in quality between the 18S and mitochondria base ML trees. Further mitochondria dataset analyses showed that the Bayesian GTR tree was more robust than the ML tree to ablation, removal of outgroups, and realignment of data subsets. The largest challenge was data processing. In particular, curating a dataset that would have better bootstrap support and posterior probabilities compared to the 18S dataset proved to be a challenge. Future directions could include testing a different clade, more and/or different genes, and other phylogenetic tree inferance software such as RAxML.</p>
</section>
<section id="references" class="level1">
<h1>References</h1>
<p>Espeland M, Breinholt J, Willmott KR, Warren AD, Vila R, Toussaint EFA, Maunsell SC, Aduse-Poku K, Talavera G, Eastwood R, Jarzyna MA, Guralnick R, Lohman DJ, Pierce NE, Kawahara AY. A Comprehensive and Dated Phylogenomic Analysis of Butterflies. <em>Curr Biol.</em> 2018 Mar 5;28(5):770-778.e5. doi: 10.1016/j.cub.2018.01.061. Epub 2018 Feb 15. PMID: 29456146.</p>
</section>

</main>
<!-- /main column -->
<script id="quarto-html-after-body" type="application/javascript">
window.document.addEventListener("DOMContentLoaded", function (event) {
  const toggleBodyColorMode = (bsSheetEl) => {
    const mode = bsSheetEl.getAttribute("data-mode");
    const bodyEl = window.document.querySelector("body");
    if (mode === "dark") {
      bodyEl.classList.add("quarto-dark");
      bodyEl.classList.remove("quarto-light");
    } else {
      bodyEl.classList.add("quarto-light");
      bodyEl.classList.remove("quarto-dark");
    }
  }
  const toggleBodyColorPrimary = () => {
    const bsSheetEl = window.document.querySelector("link#quarto-bootstrap");
    if (bsSheetEl) {
      toggleBodyColorMode(bsSheetEl);
    }
  }
  toggleBodyColorPrimary();  
  const icon = "";
  const anchorJS = new window.AnchorJS();
  anchorJS.options = {
    placement: 'right',
    icon: icon
  };
  anchorJS.add('.anchored');
  const isCodeAnnotation = (el) => {
    for (const clz of el.classList) {
      if (clz.startsWith('code-annotation-')) {                     
        return true;
      }
    }
    return false;
  }
  const onCopySuccess = function(e) {
    // button target
    const button = e.trigger;
    // don't keep focus
    button.blur();
    // flash "checked"
    button.classList.add('code-copy-button-checked');
    var currentTitle = button.getAttribute("title");
    button.setAttribute("title", "Copied!");
    let tooltip;
    if (window.bootstrap) {
      button.setAttribute("data-bs-toggle", "tooltip");
      button.setAttribute("data-bs-placement", "left");
      button.setAttribute("data-bs-title", "Copied!");
      tooltip = new bootstrap.Tooltip(button, 
        { trigger: "manual", 
          customClass: "code-copy-button-tooltip",
          offset: [0, -8]});
      tooltip.show();    
    }
    setTimeout(function() {
      if (tooltip) {
        tooltip.hide();
        button.removeAttribute("data-bs-title");
        button.removeAttribute("data-bs-toggle");
        button.removeAttribute("data-bs-placement");
      }
      button.setAttribute("title", currentTitle);
      button.classList.remove('code-copy-button-checked');
    }, 1000);
    // clear code selection
    e.clearSelection();
  }
  const getTextToCopy = function(trigger) {
      const codeEl = trigger.previousElementSibling.cloneNode(true);
      for (const childEl of codeEl.children) {
        if (isCodeAnnotation(childEl)) {
          childEl.remove();
        }
      }
      return codeEl.innerText;
  }
  const clipboard = new window.ClipboardJS('.code-copy-button:not([data-in-quarto-modal])', {
    text: getTextToCopy
  });
  clipboard.on('success', onCopySuccess);
  if (window.document.getElementById('quarto-embedded-source-code-modal')) {
    // For code content inside modals, clipBoardJS needs to be initialized with a container option
    // TODO: Check when it could be a function (https://github.com/zenorocha/clipboard.js/issues/860)
    const clipboardModal = new window.ClipboardJS('.code-copy-button[data-in-quarto-modal]', {
      text: getTextToCopy,
      container: window.document.getElementById('quarto-embedded-source-code-modal')
    });
    clipboardModal.on('success', onCopySuccess);
  }
    var localhostRegex = new RegExp(/^(?:http|https):\/\/localhost\:?[0-9]*\//);
    var mailtoRegex = new RegExp(/^mailto:/);
      var filterRegex = new RegExp('/' + window.location.host + '/');
    var isInternal = (href) => {
        return filterRegex.test(href) || localhostRegex.test(href) || mailtoRegex.test(href);
    }
    // Inspect non-navigation links and adorn them if external
 	var links = window.document.querySelectorAll('a[href]:not(.nav-link):not(.navbar-brand):not(.toc-action):not(.sidebar-link):not(.sidebar-item-toggle):not(.pagination-link):not(.no-external):not([aria-hidden]):not(.dropdown-item):not(.quarto-navigation-tool):not(.about-link)');
    for (var i=0; i<links.length; i++) {
      const link = links[i];
      if (!isInternal(link.href)) {
        // undo the damage that might have been done by quarto-nav.js in the case of
        // links that we want to consider external
        if (link.dataset.originalHref !== undefined) {
          link.href = link.dataset.originalHref;
        }
      }
    }
  function tippyHover(el, contentFn, onTriggerFn, onUntriggerFn) {
    const config = {
      allowHTML: true,
      maxWidth: 500,
      delay: 100,
      arrow: false,
      appendTo: function(el) {
          return el.parentElement;
      },
      interactive: true,
      interactiveBorder: 10,
      theme: 'quarto',
      placement: 'bottom-start',
    };
    if (contentFn) {
      config.content = contentFn;
    }
    if (onTriggerFn) {
      config.onTrigger = onTriggerFn;
    }
    if (onUntriggerFn) {
      config.onUntrigger = onUntriggerFn;
    }
    window.tippy(el, config); 
  }
  const noterefs = window.document.querySelectorAll('a[role="doc-noteref"]');
  for (var i=0; i<noterefs.length; i++) {
    const ref = noterefs[i];
    tippyHover(ref, function() {
      // use id or data attribute instead here
      let href = ref.getAttribute('data-footnote-href') || ref.getAttribute('href');
      try { href = new URL(href).hash; } catch {}
      const id = href.replace(/^#\/?/, "");
      const note = window.document.getElementById(id);
      if (note) {
        return note.innerHTML;
      } else {
        return "";
      }
    });
  }
  const xrefs = window.document.querySelectorAll('a.quarto-xref');
  const processXRef = (id, note) => {
    // Strip column container classes
    const stripColumnClz = (el) => {
      el.classList.remove("page-full", "page-columns");
      if (el.children) {
        for (const child of el.children) {
          stripColumnClz(child);
        }
      }
    }
    stripColumnClz(note)
    if (id === null || id.startsWith('sec-')) {
      // Special case sections, only their first couple elements
      const container = document.createElement("div");
      if (note.children && note.children.length > 2) {
        container.appendChild(note.children[0].cloneNode(true));
        for (let i = 1; i < note.children.length; i++) {
          const child = note.children[i];
          if (child.tagName === "P" && child.innerText === "") {
            continue;
          } else {
            container.appendChild(child.cloneNode(true));
            break;
          }
        }
        if (window.Quarto?.typesetMath) {
          window.Quarto.typesetMath(container);
        }
        return container.innerHTML
      } else {
        if (window.Quarto?.typesetMath) {
          window.Quarto.typesetMath(note);
        }
        return note.innerHTML;
      }
    } else {
      // Remove any anchor links if they are present
      const anchorLink = note.querySelector('a.anchorjs-link');
      if (anchorLink) {
        anchorLink.remove();
      }
      if (window.Quarto?.typesetMath) {
        window.Quarto.typesetMath(note);
      }
      // TODO in 1.5, we should make sure this works without a callout special case
      if (note.classList.contains("callout")) {
        return note.outerHTML;
      } else {
        return note.innerHTML;
      }
    }
  }
  for (var i=0; i<xrefs.length; i++) {
    const xref = xrefs[i];
    tippyHover(xref, undefined, function(instance) {
      instance.disable();
      let url = xref.getAttribute('href');
      let hash = undefined; 
      if (url.startsWith('#')) {
        hash = url;
      } else {
        try { hash = new URL(url).hash; } catch {}
      }
      if (hash) {
        const id = hash.replace(/^#\/?/, "");
        const note = window.document.getElementById(id);
        if (note !== null) {
          try {
            const html = processXRef(id, note.cloneNode(true));
            instance.setContent(html);
          } finally {
            instance.enable();
            instance.show();
          }
        } else {
          // See if we can fetch this
          fetch(url.split('#')[0])
          .then(res => res.text())
          .then(html => {
            const parser = new DOMParser();
            const htmlDoc = parser.parseFromString(html, "text/html");
            const note = htmlDoc.getElementById(id);
            if (note !== null) {
              const html = processXRef(id, note);
              instance.setContent(html);
            } 
          }).finally(() => {
            instance.enable();
            instance.show();
          });
        }
      } else {
        // See if we can fetch a full url (with no hash to target)
        // This is a special case and we should probably do some content thinning / targeting
        fetch(url)
        .then(res => res.text())
        .then(html => {
          const parser = new DOMParser();
          const htmlDoc = parser.parseFromString(html, "text/html");
          const note = htmlDoc.querySelector('main.content');
          if (note !== null) {
            // This should only happen for chapter cross references
            // (since there is no id in the URL)
            // remove the first header
            if (note.children.length > 0 && note.children[0].tagName === "HEADER") {
              note.children[0].remove();
            }
            const html = processXRef(null, note);
            instance.setContent(html);
          } 
        }).finally(() => {
          instance.enable();
          instance.show();
        });
      }
    }, function(instance) {
    });
  }
      let selectedAnnoteEl;
      const selectorForAnnotation = ( cell, annotation) => {
        let cellAttr = 'data-code-cell="' + cell + '"';
        let lineAttr = 'data-code-annotation="' +  annotation + '"';
        const selector = 'span[' + cellAttr + '][' + lineAttr + ']';
        return selector;
      }
      const selectCodeLines = (annoteEl) => {
        const doc = window.document;
        const targetCell = annoteEl.getAttribute("data-target-cell");
        const targetAnnotation = annoteEl.getAttribute("data-target-annotation");
        const annoteSpan = window.document.querySelector(selectorForAnnotation(targetCell, targetAnnotation));
        const lines = annoteSpan.getAttribute("data-code-lines").split(",");
        const lineIds = lines.map((line) => {
          return targetCell + "-" + line;
        })
        let top = null;
        let height = null;
        let parent = null;
        if (lineIds.length > 0) {
            //compute the position of the single el (top and bottom and make a div)
            const el = window.document.getElementById(lineIds[0]);
            top = el.offsetTop;
            height = el.offsetHeight;
            parent = el.parentElement.parentElement;
          if (lineIds.length > 1) {
            const lastEl = window.document.getElementById(lineIds[lineIds.length - 1]);
            const bottom = lastEl.offsetTop + lastEl.offsetHeight;
            height = bottom - top;
          }
          if (top !== null && height !== null && parent !== null) {
            // cook up a div (if necessary) and position it 
            let div = window.document.getElementById("code-annotation-line-highlight");
            if (div === null) {
              div = window.document.createElement("div");
              div.setAttribute("id", "code-annotation-line-highlight");
              div.style.position = 'absolute';
              parent.appendChild(div);
            }
            div.style.top = top - 2 + "px";
            div.style.height = height + 4 + "px";
            div.style.left = 0;
            let gutterDiv = window.document.getElementById("code-annotation-line-highlight-gutter");
            if (gutterDiv === null) {
              gutterDiv = window.document.createElement("div");
              gutterDiv.setAttribute("id", "code-annotation-line-highlight-gutter");
              gutterDiv.style.position = 'absolute';
              const codeCell = window.document.getElementById(targetCell);
              const gutter = codeCell.querySelector('.code-annotation-gutter');
              gutter.appendChild(gutterDiv);
            }
            gutterDiv.style.top = top - 2 + "px";
            gutterDiv.style.height = height + 4 + "px";
          }
          selectedAnnoteEl = annoteEl;
        }
      };
      const unselectCodeLines = () => {
        const elementsIds = ["code-annotation-line-highlight", "code-annotation-line-highlight-gutter"];
        elementsIds.forEach((elId) => {
          const div = window.document.getElementById(elId);
          if (div) {
            div.remove();
          }
        });
        selectedAnnoteEl = undefined;
      };
        // Handle positioning of the toggle
    window.addEventListener(
      "resize",
      throttle(() => {
        elRect = undefined;
        if (selectedAnnoteEl) {
          selectCodeLines(selectedAnnoteEl);
        }
      }, 10)
    );
    function throttle(fn, ms) {
    let throttle = false;
    let timer;
      return (...args) => {
        if(!throttle) { // first call gets through
            fn.apply(this, args);
            throttle = true;
        } else { // all the others get throttled
            if(timer) clearTimeout(timer); // cancel #2
            timer = setTimeout(() => {
              fn.apply(this, args);
              timer = throttle = false;
            }, ms);
        }
      };
    }
      // Attach click handler to the DT
      const annoteDls = window.document.querySelectorAll('dt[data-target-cell]');
      for (const annoteDlNode of annoteDls) {
        annoteDlNode.addEventListener('click', (event) => {
          const clickedEl = event.target;
          if (clickedEl !== selectedAnnoteEl) {
            unselectCodeLines();
            const activeEl = window.document.querySelector('dt[data-target-cell].code-annotation-active');
            if (activeEl) {
              activeEl.classList.remove('code-annotation-active');
            }
            selectCodeLines(clickedEl);
            clickedEl.classList.add('code-annotation-active');
          } else {
            // Unselect the line
            unselectCodeLines();
            clickedEl.classList.remove('code-annotation-active');
          }
        });
      }
  const findCites = (el) => {
    const parentEl = el.parentElement;
    if (parentEl) {
      const cites = parentEl.dataset.cites;
      if (cites) {
        return {
          el,
          cites: cites.split(' ')
        };
      } else {
        return findCites(el.parentElement)
      }
    } else {
      return undefined;
    }
  };
  var bibliorefs = window.document.querySelectorAll('a[role="doc-biblioref"]');
  for (var i=0; i<bibliorefs.length; i++) {
    const ref = bibliorefs[i];
    const citeInfo = findCites(ref);
    if (citeInfo) {
      tippyHover(citeInfo.el, function() {
        var popup = window.document.createElement('div');
        citeInfo.cites.forEach(function(cite) {
          var citeDiv = window.document.createElement('div');
          citeDiv.classList.add('hanging-indent');
          citeDiv.classList.add('csl-entry');
          var biblioDiv = window.document.getElementById('ref-' + cite);
          if (biblioDiv) {
            citeDiv.innerHTML = biblioDiv.innerHTML;
          }
          popup.appendChild(citeDiv);
        });
        return popup.innerHTML;
      });
    }
  }
});
</script>
</div> <!-- /content -->




</body></html>
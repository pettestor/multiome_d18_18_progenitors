the_plan <-
  drake_plan(

    progenitors.multiome = preprocess_multiome(),
    analysis = analyze_multiome(progenitors.multiome),
    pando = run_pando(progenitors.multiome),
    
)

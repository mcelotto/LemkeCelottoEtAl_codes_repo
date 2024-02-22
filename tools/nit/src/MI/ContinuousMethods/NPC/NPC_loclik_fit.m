
% #' Custom function for local likelihood fitting

function lfit=NPC_loclik_fit(bw,data,Grid)


[Ker_grid_Naive]=NPC_DenseNaive(bw,(data.X)',(Grid.X)');
[Ker_grid]=NPC_Kern_LL(Ker_grid_Naive,bw);

lfit.grid=Grid.X;
lfit.Kergrid=Ker_grid;
lfit.KergridNaive=Ker_grid_Naive;

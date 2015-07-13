11/03/2014 - Recuparação do solver

Este solver foi modificado com erro e recuperado para usar source-term 
e escrever um arquivo "coupling" com informaçõe da temperatura 
e densidade do fluido.

Foi feita uma simulação (myChtMultiRegionHeater) e *aparentemente* 
funciona como antes da modificação.

Este arquivo é um backup antes da alteração para que o solver seja 
capaz de manipular mais de um sólido.

12/03/2014 - Mudança - múltiplos sólidos.

Implementada a leitura do arquivo Q em todos os diretórios de sólidos.
Simulação roda.

TODO: Implementar para que sejam opcionais.
      Testar o arquivo de saída T, rho para o acoplamento.

25/03/2014 - Mudança

A leitura de termos-fonte via fvOptions foi removida no 
createSolidFields.H. O único termo-fonte desse solver é o 
implementado no solver. A leitura do arquivo ainda está 
pendente.

A remoção do fvOptions foi feita no createSolidFields.H, 
setRegionSolidFiedls e solveSolid.H/

Compilou e rodou. Necessários testes no resultado.

11/04/2014

A remoção do fvOptions (aparentemente) não afetou o 
solver.

Adicionada uma chamada a uma função extern do FORTRAN 
no arquiovo chtMultiReigionSimpleFoam.C. A função está 
implementada como uma subroutine e compilada com:

gfortran -c q.f90 [void heatflux(* float)]

o objeto q.o foi copiado no diretório 
$(WM_USER_LIBBIN) e o arquivo "options" 
dentro do Make do solver Simple 
foi modificado para inclur o objeto fortran.

A chamada funciona corretamente.

24/05/2014

Mudanças no solver. Removidos arquivos do chtMultiRegionFoam.
Gerado novo solver newChtMultiRegionSimpleFoam.

Leitura de Q externo possível de 2 formas:

1) Lendo IOobject com Q, boundaryConditions e dimensionedFields.
   Mais limpo, mas exige que o arquivo a ser lido venha no formato.

2) Alternativa: ler arquivo via fstream armazenando dados em um 
   std::vector. Dados no std::vector são depois copiados em um
   scalarField e então passados ao volScalarField.internalField() = Sf.

Ambos funcionam. Ainda não foi completamente testado com um script 
externo.

06/06/2014

Solver testado com script externo sequencialmente. 
Le e escreve Q, rho e T. Os Q lidos são efetivamente usados no cáculo.

Funcionamento sequecial: Ok.

Execução em paralelo ainda em desenvolvimento. Situação atual:
- o processo principal consegue ler o arquivo Q e recuperar os 
valores de Q internos nos processos.
- o processo principal também consegue ler os dados obtidos 
em procInfo (labelIOList) dos processos com os índices das células 
de cada domínio separado.

09/06/2014

Obtenção de Q via arquivo e cálculo com Q finalizados.
Tanto para modo sequencial quanto paralelo.
O arquivo Q no diretório 10/fuel depois do reconstructPar 
é idêntico nos dois casos.

15/04/2015 -----------------------------------------------------------

NOVO: thesisChtMultiRegionFoam

Continuação deste solver. A ser usado na tese.
Mapping 1 para 1.

Testado inicialmente na buril-lx.

15/04/2015

Leitura dos arquivos em constant/set para uso no mapeamento dos campos
em cada processador para um vetor global com dados de rho, T e Q.

30/04/2015

Os sets relacionando as células na região ao mesh completo são lidos
via IOobjects.
Logo após, uma lista de labelList é preenchida com os valores da hashtable
original cellSets.

Estas listas contém, na posição atual da malha, qual a célula correspondente
na malha completa. Isso feito por região.

TODO: Ler valores de T, rho e Q para o vetor geral
      Comparar as posições com o que o paraview vai dar.
      
10/05/2015

Os sets são modificados via chamada de sistema sed para que a classe
no arquivo 'cellSet' seja substituída por 'labelList'. Isso é necessário
pois cellSets são automaticamente lidos em uma hashtable e os dados
perdem o número relativo de posição. Com a mudança, passa a existir uma
lista com as posições de cada célula na malha.

Leitura de valores:

Sequencial: lê por região e escreve na posição relativa na malha completa.
Paralelo: lê por processador, escreve a posição relativa da célula em cada
	  processador no vetor da região. Depois, escreve a posição relativa
	  por região no vetor da malha completa.

TODO: Colocar valores de Q nos processadores.
      Normalizar neutrônica.
      Modificar neutrônica.

09/07/2015

Sequencial e paralelo implementados.
A neutronica é testada por um vetor com valores aleatórios entre 0.9 e 1.1.
Este valor é normalizado no termo-fonte.

A neutronica é executada a cada 100 iterações da termohidráulica. Está fixo.

TODO: Criar IOobject para controlar as iteracoes.
      Adaptar a neutronica verdadeira e acrescentar.
      

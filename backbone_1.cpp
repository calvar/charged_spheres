//Paralelismo grueso MEJOR

// N = 100000




//Arreglo posicion 3*N*double
//Arreglo velocidad 3*N*double
//Arreglo fuerza 3*N*double
//Arreglo de masa N*double
//Arreglo de carga N*double
//Algunas variables para acumular

//11*N*8 * Bytes = 1x10**7 Bytes

//PCI 8 Gtrans / s

//P100 -> 16GB
//K80 -> 16GB
//M40 -> 24GB
//K20 -> 5GB


//Copiar memoria


///Parte REAL--------------------------------
for(int i = 0; i < N; ++i){

  //Version mejor, segun yo mero
  int j = i+1 - N*static_cast<int>(floor((i+1)/N + 0.5));
  int cnt = 0;
  while(cnt < mx){
    //Mirar que operaciones se pueden vectorizar dentro de este ciclo

    //Hay if --> problema para la vectorizacion!!!!!!
    
    print("hola");
    //wdwdwdw

    //Que va aqui:
    //10-20 operaciones basicas ( + - * / )
    //un par de funciones especiales ()
    
  }
  j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
  cnt++;

}

//Copiar de vuelta


///Parte reciproca--------------------------



for(int kn = 0; kn < K_max3; ++kn){ //K_max_3 es del orden de 9**3
  //Hay if

  //Copiar memoria
  
  for(int i = 0; i < N; ++i){  //N es del orden de 10000+
    //varias operaciones de punto flotante
  }

  //varias operaciones entre elementos de un arreglo
  
  for(int i = 0; i < N; ++i){
    //Varias operaciones de punto flotante
  }

  //copiar memoria
}





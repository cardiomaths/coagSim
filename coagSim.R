coagSim <- function(donors, inhib) {
  
  if (donors != "All" && donors != "Female" && donors != "Male") {
    stop("The first argument must be All, Male or Female")
    }

  if (inhib != 'none' && inhib != 'Xa' && inhib != 'Th' && inhib != 'Warfarin25') {
    stop("The second argument must be none, Xa, Th or Warfarin25")
  }
  
    library(deSolve)
  
  times <-
    seq(from = 0, to = 10000, by = 0.02)# minimum tolerance to work
  
  XaI0 <- 0
  ThI0 <- 0
  
  if (inhib == "Xa") {
    XaI0 <- 1
    ThI0 <- 0
  }
  
  if (inhib == "Th") {
    XaI0 <- 0
    ThI0 <- 1
  }
  
  parameters <- c(
    a1 <- 3.1e-3,
    a2 <- 3.2e6,
    a3 <- 3.1e-3,
    a4 <- 2.3e7,
    a5 <- 4.4e5,
    a6 <- 1.3e7,
    a7 <- 2.3e4,
    a8 <- 1.05,
    a9 <- 2.5e7,
    a10 <- 6,
    a11 <- 19,
    a12 <- 2.2e7,
    a13 <- 2.4,
    a14 <- 1e7,
    a15 <- 1.8,
    a16 <- 7.5e3,
    a17 <- 2e7,
    a18 <- 5e-3,
    a19 <- 1e7,
    a20 <- 1e-3,
    a21 <- 1e8,
    a22 <- 8.2,
    a23 <- 2.2e4,
    a24 <- 6e-3,
    a25 <- 1e-3,
    a26 <- 2e7,
    a27 <- 0.2,
    a28 <- 4e8,
    a29 <- 103,
    a30 <- 1e8,
    a31 <- 63.5,
    a32 <- 2.3e8,
    #altered in later papers
    a33 <- 3.6e-4,
    a34 <- 9e5,
    a35 <- 1.1e-4,
    a36 <- 3.2e8,
    a37 <- 5e7,
    a38 <- 4.2e3,
    #altered in later papers
    a39 <- 7.1e3,
    a40 <- 4.9e2,
    a41 <- 7.1e3,
    a42 <- 2.3e2,
    a43 <- 5.7e3,
    #newly introduced in later papers
    a44 <- 3e6,
    XaI <- 0,
    ThI <- 0
  )
  
  u3 = 1.6e-7 #x 160nM
  u4 = 2e-8 #v 20nM
  u5 = 1.4e-6 #Thrombin 1400nM
  u6 = 3.4e-6 #ATIII 3400nM
  u8 = 9e-8 #ix 90nM
  u9 = 7e-10 #viii 0.7nM
  u32 = 1.39e-9 #TFPI 2.4nM
  u36 = 1e-8 #VII 10nM
  u35 = 1e-10 #VIIa 0.1nM
  u34 = 1e-12
  
  if (inhib == 'Warfarin25') {
    u3 <- u3 * 0.25
    u5 <- u5 * 0.25
    u8 <- u8 * 0.25
    u36 <- u36 * 0.25
    u35 <- u35 * 0.25
  }
  
  ini <- c(
    Z1 = 0,
    Z2 = 0,
    Z3 = u3,
    Z4 = u4,
    Z5 = u5,
    Z6 = u6,
    Z7 = 0,
    Z8 = u8,
    Z9 = u9,
    Z10 = 0,
    Z11 = 0,
    Z12 = 0,
    Z13 = 0,
    Z14 = 0,
    Z15 = 0,
    Z16 = 0,
    Z17 = 0,
    Z18 = 0,
    Z19 = 0,
    Z20 = 0,
    Z21 = 0,
    Z22 = 0,
    Z23 = 0,
    Z24 = 0,
    Z25 = 0,
    Z26 = 0,
    Z27 = 0,
    Z28 = 0,
    Z29 = 0,
    Z30 = 0,
    Z31 = 0,
    Z32 = u32,
    Z33 = 0,
    Z34 = u34,
    Z35 = u35,
    Z36 = u36,
    Z37 = 0,
    Z38 = 0,
    Z39 = 0,
    Z40 = 0,
    Z41 = 0,
    Z42 = 0,
    Z43 = 0
  )
  
  modelH <- function (t, ini, parmameters) {
    with(as.list(c(ini, parameters)), {
      dZ34 <- -a2 * Z34 * Z36 + a1 * Z37 - a4 * Z34 * Z35 + a3 * Z12 #TF
      dZ36 <-
        -a2 * Z34 * Z36 + a1 * Z37 - a5 * Z12 * Z36 - a6 * Z13 * Z36 - a7 * Z16 *
        Z36 #VII
      dZ37 <- -a1 * Z37 + a2 * Z34 * Z36 #TF:VII
      dZ35 <-
        -a4 * Z34 * Z35 + a3 * Z12 + a5 * Z12 * Z36 + a6 * Z13 * Z36 + a7 * Z16 *
        Z36 #VIIa
      dZ12 <-
        -a3 * Z12 + a4 * Z34 * Z35 - a9 * Z12 * Z3 + a8 * Z28 - a12 * Z12 * Z13 +
        a11 * Z29 - a14 * Z12 * Z8 + a13 * Z38 + a15 * Z38 - a37 * Z12 * Z31 - a42 *
        Z12 * Z6 #TF:VIIa Extrinsic tenase (n)
      dZ28 <- a9 * Z12 * Z3 - a10 * Z28 - a8 * Z28   #TF:VIIa:X
      dZ29 <-
        a10 * Z28 + a12 * Z12 * Z13 - a11 * Z29 - a36 * Z29 * Z32 + a35 * Z33 #TF:VIIa:Xa
      dZ38 <- a14 * Z12 * Z8 - a13 * Z38 - a15 * Z38  #TF:VIIa:IX
      
      dZ32 <- -a34 * Z13 * Z32 + a33 * Z31 - a36 * Z29 * Z32 + a35 * Z33 #TFPI
      dZ31 <- a34 * Z13 * Z32 - a33 * Z31 - a37 * Z12 * Z31 #b_xa p
      dZ33 <- a36 * Z29 * Z32 - a35 * Z33 + a37 * Z12 * Z31 #VIIa:TF:Xa:TFPI
      
      dZ6  <-
        -a38 * Z13 * Z6 - a39 * Z17 * Z6 - a40 * Z18 * Z6 - a41 * Z16 * Z6 - a42 *
        Z12 * Z6 #AT
      dZ22 <- a38 * Z13 * Z6 #Xa:AT
      dZ25 <- a39 * Z17 * Z6 #mIIa:AT
      dZ26 <- a40 * Z18 * Z6 #IXa:AT
      dZ24 <- a41 * Z16 * Z6 #IIa:AT
      dZ7  <- a42 * Z12 * Z6 #TF:VIIa:AT
      
      dZ3  <- -a9 * Z12 * Z3 + a8 * Z28 - a21 * Z20 * Z3 + a20 * Z27 + a25 *
        Z27 #X
      dZ4  <- -a26 * Z16 * Z4 - a44 * Z17 * Z4 #V
      dZ5  <- -a16 * Z13 * Z5 - a30 * Z15 * Z5 + a29 * Z10 #prothrombin
      dZ13 <-
        -a12 * Z12 * Z13 + a11 * Z29 - a28 * Z13 * Z14 + a27 * Z15 + a22 * Z27 -
        a34 * Z13 * Z32 + a33 * Z31 - a38 * Z13 * Z6  + a43 * Z18 * Z3 - XaI * 1e8 *
        Z13 * Z11 + XaI * 0.04 * Z30 #Xa
      dZ14 <- a26 * Z16 * Z4 - a28 * Z13 * Z14 + a27 * Z15  + a44 * Z17 * Z4 #Va
      dZ15 <-
        a28 * Z13 * Z14 - a27 * Z15 - a30 * Z15 * Z5 + a29 * Z10 + a31 * Z10 - XaI *
        1e8 * Z15 * Z11 + XaI * 0.21 * Z1 #Prothrombinase
      dZ10 <- a30 * Z15 * Z5 - a29 * Z10 - a31 * Z10 #Xa:Va:II
      dZ16 <-
        a16 * Z13 * Z5 + a32 * Z17 * Z15 - a41 * Z16 * Z6  - ThI * 2.5 * 1e6 *
        Z41 * Z16 + 0.05 * Z42 #Thrombin
      dZ17 <-
        a31 * Z10 - a32 * Z17 * Z15 - a39 * Z17 * Z6  - ThI * 8.5 * 1e4 * Z41 *
        Z17 - 0.05 * Z43 #mIIa
      
      dZ8  <- -a14 * Z12 * Z8 + a13 * Z38 #IX
      dZ9  <- -a17 * Z16 * Z9 #VIII
      dZ18 <-
        a15 * Z38 - a19 * Z19 * Z18 + a18 * Z20 + a25 * Z27 + a25 * Z20 - a40 *
        Z18 * Z6 #IXa
      dZ19 <-
        a17 * Z16 * Z9 - a19 * Z19 * Z18 + a18 * Z20 - a24 * Z19 + a23 * Z21 *
        Z23 #VIIIa
      dZ21 <- a24 * Z19 + a25 * Z27 + a25 * Z20 - a23 * Z21 * Z23 #VIIIa1
      dZ23 <- a24 * Z19 + a25 * Z27 + a25 * Z20 - a23 * Z21 * Z23 #VIIIa2
      
      dZ20 <-
        a19 * Z19 * Z18 - a18 * Z20 - a21 * Z20 * Z3 + a20 * Z27 + a22 * Z27 - a25 *
        Z20 #Intrinsic tenase
      dZ27 <-
        a21 * Z20 * Z3 - a20 * Z27 - a22 * Z27 - a25 * Z27 #intrinsic tenase bound to X
      
      dZ39 <- Z16 + Z17
      dZ40 <- Z13
      
      dZ11 <-
        -XaI * 1e8 * Z13 * Z11 + XaI * 0.04 * Z30 - XaI * 1e8 * Z15 * Z11 + XaI *
        0.21 * Z1
      dZ30 <- XaI * 1e8 * Z13 * Z11 - XaI * 0.04 * Z30
      dZ1 <- XaI * 1e8 * Z15 * Z11 - XaI * 0.21 * Z1
      dZ2 <- 0
      
      
      #Factor IIa inhibitor (DAPA)
      dZ41 <-
        -ThI * 2.5 * 1e6 * Z41 * Z16 + 0.05 * Z42 - ThI * 8.5 * 1e4 * Z41 * Z17 + 0.05 *
        Z43
      dZ42 <- ThI * 2.5 * 1e6 * Z41 * Z16 - 0.05 * Z42
      dZ43 <- ThI * 8.5 * 1e4 * Z41 * Z17 - 0.05 * Z43
      
      # return
      list(
        c(
          dZ1,
          dZ2,
          dZ3,
          dZ4,
          dZ5,
          dZ6,
          dZ7,
          dZ8,
          dZ9,
          dZ10,
          dZ11,
          dZ12,
          dZ13,
          dZ14,
          dZ15,
          dZ16,
          dZ17,
          dZ18,
          dZ19,
          dZ20,
          dZ21,
          dZ22,
          dZ23,
          dZ24,
          dZ25,
          dZ26,
          dZ27,
          dZ28,
          dZ29,
          dZ30,
          dZ31,
          dZ32,
          dZ33,
          dZ34,
          dZ35,
          dZ36,
          dZ37,
          dZ38,
          dZ39,
          dZ40,
          dZ41,
          dZ42,
          dZ43
        )
      )
    })
  }
  
  df  <- read.table("coagData.csv", sep = ",", header = T)
  
  if (donors == "All") {
    caTF0 <-   df[1, 2]
    caTFPI0 <- df[1, 3]
    caII0 <-   df[1, 4] / 100
    caV0 <-    df[1, 5] / 100
    caVII0 <-  df[1, 6] / 100
    caX0 <-    df[1, 7] / 100
    caVIII0 <- df[1, 8] / 100
    caIX0 <-   df[1, 9] / 100
    caAT0 <-   df[1, 10] / 100
    
    #Load Control means
    coTF0 <-   df[2, 2]
    coTFPI0 <- df[2, 3]
    coII0 <-   df[2, 4] / 100
    coV0 <-    df[2, 5] / 100
    coVII0 <-  df[2, 6] / 100
    coX0 <-    df[2, 7] / 100
    coVIII0 <- df[2, 8] / 100
    coIX0 <-   df[2, 9] / 100
    coAT0 <-   df[2, 10] / 100
  }
  
  if (donors == "Male") {
    caTF0 <-   df[3, 2]
    caTFPI0 <- df[3, 3]
    caII0 <-   df[3, 4] / 100
    caV0 <-    df[3, 5] / 100
    caVII0 <-  df[3, 6] / 100
    caX0 <-    df[3, 7] / 100
    caVIII0 <- df[3, 8] / 100
    caIX0 <-   df[3, 9] / 100
    caAT0 <-   df[3, 10] / 100
    
    #Load Control means
    coTF0 <-   df[4, 2]
    coTFPI0 <- df[4, 3]
    coII0 <-   df[4, 4] / 100
    coV0 <-    df[4, 5] / 100
    coVII0 <-  df[4, 6] / 100
    coX0 <-    df[4, 7] / 100
    coVIII0 <- df[4, 8] / 100
    coIX0 <-   df[4, 9] / 100
    coAT0 <-   df[4, 10] / 100
  }
  
  if (donors == "Female") {
    caTF0 <-   df[5, 2]
    caTFPI0 <- df[5, 3]
    caII0 <-   df[5, 4] / 100
    caV0 <-    df[5, 5] / 100
    caVII0 <-  df[5, 6] / 100
    caX0 <-    df[5, 7] / 100
    caVIII0 <- df[5, 8] / 100
    caIX0 <-   df[5, 9] / 100
    caAT0 <-   df[5, 10] / 100
    
    #Load Control means
    coTF0 <-   df[6, 2]
    coTFPI0 <- df[6, 3]
    coII0 <-   df[6, 4] / 100
    coV0 <-    df[6, 5] / 100
    coVII0 <-  df[6, 6] / 100
    coX0 <-    df[6, 7] / 100
    coVIII0 <- df[6, 8] / 100
    coIX0 <-   df[6, 9] / 100
    coAT0 <-   df[6, 10] / 100
  }
  
  iniCa <- c(
    Z1 = 0,
    Z2 = 0,
    Z3 = u3 * caX0,
    #x
    Z4 = u4 * caV0,
    #v
    Z5 = u5 * caII0,
    #Prothrombin
    Z6 = u6 * caAT0,
    #ATIII
    Z7 = 0,
    Z8 = u8 * caIX0,
    #ix
    Z9 = u9 * caVIII0,
    #viii
    Z10 = 0,
    Z11 = 0,
    Z12 = 0,
    Z13 = 0,
    Z14 = 0,
    Z15 = 0,
    Z16 = 0,
    Z17 = 0,
    Z18 = 0,
    Z19 = 0,
    Z20 = 0,
    Z21 = 0,
    Z22 = 0,
    Z23 = 0,
    Z24 = 0,
    Z25 = 0,
    Z26 = 0,
    Z27 = 0,
    Z28 = 0,
    Z29 = 0,
    Z30 = 0,
    Z31 = 0,
    Z32 = u32 * caTFPI0,
    Z33 = 0,
    Z34 = u34 * caTF0,
    #TF
    Z35 = u35 * caVII0,
    Z36 = u36 * caVII0,
    #%VII
    Z37 = 0,
    Z38 = 0,
    Z39 = 0,
    Z40 = 0,
    Z41 = 0,
    Z42 = 0,
    Z43 = 0
  )
  
  outCa <-
    ode(
      y = iniCa,
      times = times,
      func = modelH,
      parms = parameters
    )
  
  iniCo <- c(
    Z1 = 0,
    Z2 = 0,
    Z3 = 1.6e-7 * coX0,
    #x
    Z4 = 2e-8 * coV0,
    #v
    Z5 = 1.4e-6 * coII0,
    #Prothrombin
    Z6 = 3.4e-6 * coAT0,
    #ATIII
    Z7 = 0,
    Z8 = 9e-8 * coIX0,
    #ix
    Z9 = 7e-10 * coVIII0,
    #viii
    Z10 = 0,
    Z11 = 0,
    Z12 = 0,
    Z13 = 0,
    Z14 = 0,
    Z15 = 0,
    Z16 = 0,
    Z17 = 0,
    Z18 = 0,
    Z19 = 0,
    Z20 = 0,
    Z21 = 0,
    Z22 = 0,
    Z23 = 0,
    Z24 = 0,
    Z25 = 0,
    Z26 = 0,
    Z27 = 0,
    Z28 = 0,
    Z29 = 0,
    Z30 = 0,
    Z31 = 0,
    Z32 = u32 * coTFPI0,
    #TFPI
    Z33 = 0,
    Z34 = u34 * coTF0,
    #TF
    Z35 = u35 * coVII0,
    #VIIa
    Z36 = u36 * coVII0,
    #%VII
    Z37 = 0,
    Z38 = 0,
    Z39 = 0,
    Z40 = 0,
    Z41 = 0,
    Z42 = 0,
    Z43 = 0
  )
  
  outCo <-
    ode(
      y = iniCo,
      times = times,
      func = modelH,
      parms = parameters
    )
  
  png(
    filename = paste("Fig_ProcoagR",donors,inhib,".png",sep=""),
    width = 12,
    height = 8,
    units = 'in',
    res = 400
  )
  
  par(mfrow = c(3, 3))
  par(oma = c(0.5, 0, 0.5, 0.5)) #Sets outer region around plot (bottom, left, top, right)
  par(mgp = c(1.1, 0.5, 0)) #Sets space between data area and axis labels
  par(mar = c(3, 2, 1, 1.75)) #sets the space between plots
  par(
    ps = 8,
    cex = 1,
    cex.main = 1,
    cex.axis = 0.8,
    cex.lab = 0.8,
    font.main = 1
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 35]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "TF (nM)",
       xlim = c(0.0, 20),
       ylim = c(0, 0.015),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 35]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 6]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "II (nM)",
       xlim = c(0.0, 20),
       ylim = c(0, 1800),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 6]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 5]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "V (nM)",
       xlim = c(0.0, 20),
       ylim = c(0, 28),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 5]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 37]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "VII (nM)",
       xlim = c(0.0, 60),
       ylim = c(0, 13),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 37]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 10]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "VIII (nM)",
       xlim = c(0.0, 60),
       ylim = c(0, 1.5),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 10]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 9]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "IX (nM)",
       xlim = c(0.0, 100),
       ylim = c(0, 150),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 9]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 4]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "X (nM)",
       xlim = c(0.0, 100),
       ylim = c(0, 220),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 4]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 33]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "TFPI (nM)",
       xlim = c(0.0, 100),
       ylim = c(0, 1.50),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 33]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 7]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "AT (nM)",
       xlim = c(0.0, 100),
       ylim = c(0, 4000),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 7]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  
  dev.off()
  
  
  png(
    paste("Fig_ActiveR",donors,inhib,".png",sep=""),
    width = 12,
    height = 8,
    units = 'in',
    res = 400
  )
  
  par(mfrow = c(4, 4))
  par(oma = c(0.5, 0, 0.5, 0.5)) #Sets outer region around plot (bottom, left, top, right)
  par(mgp = c(1.1, 0.5, 0)) #Sets space between data area and axis labels
  par(mar = c(3, 2, 1, 1.75)) #sets the space between plots
  par(
    ps = 8,
    cex = 1,
    cex.main = 1,
    cex.axis = 0.8,
    cex.lab = 0.8,
    font.main = 1
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 17] + outCa[, 18]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Thrombin (nM)",
       xlim = c(0.0, 20.0),
       ylim = c(0, 750),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 17] + outCo[, 18]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 17] + outCa[, 18]) / max(outCo[, 17] + outCo[, 18]) -
                     1) * 100, digits = 2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 17] + outCa[, 18]) / which.max(outCo[, 17] + outCo[, 18]) -
        1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 17]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Thrombin (nM)",
       xlim = c(0.0, 20.0),
       ylim = c(0, 750),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 17]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 17]) / max(outCo[, 17]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 17]) / which.max(outCo[, 17]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 15]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor Va (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 25),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 15]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 15]) / max(outCo[, 15]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 15]) / which.max(outCo[, 15]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    bty = "n",
    box.lwd = 3,
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 36]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor VIIa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 25),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 36]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 36]) / max(outCo[, 36]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 36]) / which.max(outCo[, 36]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    bty = "n",
    box.lwd = 3,
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 20]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor VIIIa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 2),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 20]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 20]) / max(outCo[, 20]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 20]) / which.max(outCo[, 20]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 22] + outCa[, 24]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor VIIIa, other (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 3),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 22] + outCo[, 24]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 22] + outCa[, 24]) / max(outCo[, 22] + outCo[, 24]) -
                     1) * 100, digits = 2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 22] + outCa[, 24]) / which.max(outCo[, 22] + outCo[, 24]) -
        1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 19]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor IXa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.10),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 19]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 19]) / max(outCo[, 19]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 19]) / which.max(outCo[, 19]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 14]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor Xa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 20),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 14]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 14]) / max(outCo[, 14]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 14]) / which.max(outCo[, 14]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 13]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor TF:VIIa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.0003),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 13]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 13]) / max(outCo[, 13]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 13]) / which.max(outCo[, 13]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 29]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor TF:VIIa:X (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.0003),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 29]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 29]) / max(outCo[, 29]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 29]) / which.max(outCo[, 29]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 39]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor TF:VIIa:IX (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.0003),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 39]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 39]) / max(outCo[, 39]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 39]) / which.max(outCo[, 39]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 16] + outCa[, 11]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor Xa:Va (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 25),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 16] + outCo[, 11]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 16] + outCa[, 11]) / max(outCo[, 16] + outCo[, 11]) -
                     1) * 100, digits = 2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 16] + outCa[, 11]) / which.max(outCo[, 16] + outCo[, 11]) -
        1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 21] + outCa[, 28]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "Factor IXa:VIIIa (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.01),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 21] + outCo[, 28]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 21] + outCa[, 28]) / max(outCo[, 21] + outCo[, 28]) -
                     1) * 100, digits = 2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 21] + outCa[, 28]) / which.max(outCo[, 21] + outCo[, 28]) -
        1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 23] + outCa[, 26] + outCa[, 27] + outCa[, 25] + outCa[, 28]) *
         1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "AT usage (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 2000),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 23] + outCo[, 26] + outCo[, 27] + outCo[, 25] + outCo[, 28]) *
          1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((
      max(outCa[, 23] + outCa[, 26] + outCa[, 27] + outCa[, 25] + outCa[, 28]) /
        max(outCo[, 23] + outCo[, 26] + outCo[, 27] + outCo[, 25] + outCo[, 28]) -
        1
    ) * 100, digits = 2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 23] + outCa[, 26] + outCa[, 27] + outCa[, 25] + outCa[, 28]) /
        which.max(outCo[, 23] + outCo[, 26] + outCo[, 27] + outCo[, 25] + outCo[, 28]) -
        1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 32]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "TFPI usage (Xa) (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 1.5),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 32]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 32]) / max(outCo[, 32]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 32]) / which.max(outCo[, 32]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  plot((outCa[, 1]) / 60,
       (outCa[, 34]) * 1e9,
       type = "l",
       lwd = 2,
       lty = 1,
       xlab = "time (mins)",
       ylab = "TFPI usage (full) (nM)",
       xlim = c(0.0, 100.0),
       ylim = c(0, 0.1),
       col = "darkred"
  )
  lines((outCo[, 1]) / 60,
        (outCo[, 34]) * 1e9,
        type = "l",
        lwd = 2,
        lty = 2,
        col = "darkgreen"
  )
  legend(
    "topright",
    c(paste(round((max(outCa[, 34]) / max(outCo[, 34]) - 1) * 100, digits =
                    2
    ), "% higher"), paste(round((abs(
      which.max(outCa[, 34]) / which.max(outCo[, 34]) - 1
    )) * 100, digits = 2
    ), "% earlier")),
    text.col = "black",
    bty = "n",
    title = "Case",
    cex = 0.8
  )
  
  dev.off()
  
  return("Done")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
  }

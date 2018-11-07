// Mathias May

import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;
import mm82.Math.*;
import static mm82.Math.FunctionValues.*;

import java.awt.image.BufferedImage;
//import java.io.File;
//import java.io.IOException;
import java.net.URL;
import javax.imageio.ImageIO;
import java.awt.Color;

public class MinimizeFunction
{



    private static boolean printing_out_debug_informations_is_activated = false;
    
    public static void main(String[] args) throws IOException {
        
        int n;
        int k;
        
        long TimeBefore;
        long TimeAfter;
        long TimeSoFar;
        
        int n_X;
        int n_Y;
        
        int IterationNumber=1;
        
        long LastPrintoutTime  = System.nanoTime();
        long LastImageSaveTime = System.nanoTime();
        long LastImageSaveTimeRhos = System.nanoTime();
        long LastImageSaveTimeAbsCs = System.nanoTime();
        
        /* Parameters for plotting */
        int MarginLeft = 50;
        int MarginRight = 50;
        int MarginTop = 50;
        int MarginBottom = 50;
        
        int BoxLineThickness=3;
        
        /* */
        int NumberOfSites=3;
        int MinNumberOfParticles=0;
        int MaxNumberOfParticles=13;
        int NumberOfCoordinatesPerSite=MaxNumberOfParticles-MinNumberOfParticles+1;
        int NumberOfParameters=2*NumberOfSites*NumberOfCoordinatesPerSite;
                
                
        BufferedImage image;
        BufferedImage[] image_rho = new BufferedImage[NumberOfSites];
        BufferedImage[][] image_absC = new BufferedImage[NumberOfCoordinatesPerSite][NumberOfSites];
        
        double[] rho = new double[NumberOfSites];
        double[][] absC = new double[NumberOfSites][NumberOfCoordinatesPerSite];
        
        double[] parameters = new double[NumberOfParameters];
        double[] BestSoFar = new double[NumberOfParameters];
        
        int CounterForTotalNumberOfMinimizationsPerPixelInPhaseDiagram = 0;
        
        double ToPlotYMin = -1.0;
        double ToPlotYMax =  0.001;
        
        int NumberOfPixelAlongXAxis =  100;
        int NumberOfPixelAlongYAxis = 800;
        
        int NumberOfMinimizationsPerPixelInPhaseDiagram = 200;//10
        int NumberOfSuccessfulMinimizationsPerPixelInPhaseDiagramAfterWhichToStop = 10;// OK: 5    ; better: 10-100
        
        int NumberOfIterations = 20000;//(int)Math.pow(10.0, 5)
        int MinNumberOfIterations = 100000;//9000
        int MaxNumberOfIterations = 1000000;//100000 (if NumberOfMinimizationsPerPixelInPhaseDiagram == 100)
        double a0 =
            (
                   Math.log((double)MinNumberOfIterations)
                  *
                   (double)NumberOfMinimizationsPerPixelInPhaseDiagram
                -
                   Math.log((double)MaxNumberOfIterations)
                  *
                   1.0
            )
            /
            (
                 (double)NumberOfMinimizationsPerPixelInPhaseDiagram
                -1.0
            )
        ;
        double a1 =
            (
                 Math.log((double)MaxNumberOfIterations)
                -
                 Math.log((double)MinNumberOfIterations)
            )
            /
            (
                 (double)NumberOfMinimizationsPerPixelInPhaseDiagram
                -1.0
            )
        ;
        
        System.out.println( "a0 = " + a0 );
        System.out.println( "a1 = " + a1 );
        
        double   GoalYDifferenceToLastStep = Math.pow(10.0, -4); //Math.pow(10.0, -12)
        double[] GoalXDifferenceToLastStep = new double[NumberOfParameters];
            Arrays.fill( GoalXDifferenceToLastStep, Math.pow(10.0, -4) ); //Math.pow(10.0, -12)
        /*
        double   GoalYRelativeError = Math.pow(10.0, -6); //Math.pow(10.0, -6)
        double[] GoalXRelativeError = new double[NumberOfParameters];
            Arrays.fill( GoalXRelativeError, Math.pow(10.0, -6) ); //Math.pow(10.0, -6)
        */
        
        double V = 2.0;
        
        //double J_min  = -0.5 - 1.0/16.0 ;
        double J_min  = -1.0 ;  //-0.25
        //double J_max  =  0.0 + 1.0/16.0  ;
        double J_max  =  1.0  ;  //-0.2
        
        double J_step =  (J_max-J_min)/(double)(NumberOfPixelAlongXAxis-1);
        
        
        //double mu_min  = -0.5 - 1.0/16.0 ;
        double mu_min  =  0.0-2.0 ;  //2.22
        //double mu_max  =  0.0 + 1.0/16.0 ;
        double mu_max  =  12.0+2.0 ;  //2.24
        
        double mu_step =  (mu_max-mu_min)/(double)(NumberOfPixelAlongYAxis-1);
        
        double J_for;
        double mu_for;
        
        double temp;
        
        int red;
        int green;
        int blue;
        int alpha;
        
        float hue;
        float saturation;
        float value;
        
        int type;
        Color ColorObject;
        Energy EnergyFunction = new Energy(MinNumberOfParticles, MaxNumberOfParticles);

        
        java.util.Date date = new java.util.Date();
        
        double PercentageSoFar=0.0;
        double EstimatedTotalTime=0.0;
        double EstimatedTimeLeft=0.0;
        
        //double[] StartValuesx = new double[NumberOfParameters];
        double[] StepSizesx   = new double[NumberOfParameters];
        double[] StepSizeModificationx   = new double[NumberOfParameters];
        double   StepSizeModificationDecreaseValue   = 4.0; //4.0
        double   StepSizeModificationIncreaseValue   = 1.1; //1.1
        double   StepSizeModificationDefaultValue   = 0.5;  //0.5
        double   StepSizeModificationMinValue   = Math.pow(10.0, -20);
        double   StepSizeModificationMaxValue   = 0.5;
        
        double LowestEnergy=0.0;
        double[] Energies = new double[NumberOfMinimizationsPerPixelInPhaseDiagram];
        
        double[] gradient;
        
        boolean all_goals_are_fulfilled = false;
        boolean all_goals_are_fulfilled_involving_absolute_errors = false;
        boolean all_goals_are_fulfilled_involving_relative_errors = false;
        boolean all_goals_are_fulfilled_in_at_least_one_for_current_pixel = false;
        int NumberOfSuccessesForCurrentPixel=0;
        
        double LastEnergy = 0.0;
        double CurrentEnergy = 0.0;
        double BestSoFarInnerLoopEnergy = 0.0;
        double[] LastX    = new double[NumberOfParameters];
            Arrays.fill( LastX, 0.0 ); //0.0000001
        double[] CurrentX = new double[NumberOfParameters];
            Arrays.fill( CurrentX, 0.0 ); //0.0000001
        double[] BestSoFarInnerLoopX = new double[NumberOfParameters];
            Arrays.fill( BestSoFarInnerLoopX, 0.0 ); //0.0000001
            
        double DoubleReadFromFile = 0.0;
        int       IntReadFromFile = 0;
        
        
        DataOutputStream OutputFileStreamWithFormatedData_Current = new DataOutputStream( new FileOutputStream( "./output/current/dataset.dat" ) );
        
        DataOutputStream OutputFileStreamWithFormatedData = new DataOutputStream( new FileOutputStream( "./output/dataset " + date.toString() + ".dat" ) );
        
        OutputFileStreamWithFormatedData.writeInt( NumberOfPixelAlongXAxis );
        OutputFileStreamWithFormatedData.writeInt( NumberOfPixelAlongYAxis );
        
        OutputFileStreamWithFormatedData_Current.writeInt( NumberOfPixelAlongXAxis );
        OutputFileStreamWithFormatedData_Current.writeInt( NumberOfPixelAlongYAxis );
        
        DataInputStream InputFileStreamWithFormatedData_Current = new DataInputStream( new FileInputStream( "./output/current/dataset.dat" ) );
        
        DataInputStream InputFileStreamWithFormatedData = new DataInputStream( new FileInputStream( "./output/dataset " + date.toString() + ".dat" ) );
        
        /*
        while (true) {
            try {
                temp = InputFileStreamWithFormatedData.readDouble();
                System.out.println(temp);
            }
            catch ( java.io.EOFException e ) {
                InputFileStreamWithFormatedData.close();
                break;
            }
        }
        */
        
                if (true) {
                
        printing_out_debug_informations_is_activated = false;
        
        EnergyFunction.SetInterSiteInteractionStrength(V);

        type = BufferedImage.TYPE_4BYTE_ABGR;
        ColorObject = new Color(255, 0, 0, 0);

        image = new BufferedImage(
            NumberOfPixelAlongXAxis + MarginLeft + MarginRight,
            NumberOfPixelAlongYAxis + MarginTop  + MarginBottom,
            type
        );

        for(n_X = 0; n_X<=NumberOfPixelAlongXAxis-1; n_X++) {
            for(n_Y = 0; n_Y<=NumberOfPixelAlongYAxis-1; n_Y++) {
                ColorObject = new Color(
                    (int)Math.floor(Math.random()*(double)256),
                    (int)Math.floor(Math.random()*(double)256),
                    (int)Math.floor(Math.random()*(double)256),
                    255
                );
                image.setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
            }
        }
        
        for(n_X = -BoxLineThickness; n_X<=NumberOfPixelAlongXAxis-1+BoxLineThickness; n_X++) {
            for(n_Y = 0; n_Y<=BoxLineThickness-1; n_Y++) {
                ColorObject = new Color(
                    (int)0,
                    (int)0,
                    (int)0,
                    255
                );
                image.setRGB(n_X + MarginLeft, n_Y +MarginTop-1 -BoxLineThickness+1, ColorObject.getRGB());
                image.setRGB(n_X + MarginLeft, n_Y +NumberOfPixelAlongYAxis+1 +MarginTop-1, ColorObject.getRGB());
            }
        }
        
        for(n_X = 0; n_X<=BoxLineThickness-1; n_X++) {
            for(n_Y = -BoxLineThickness; n_Y<=NumberOfPixelAlongYAxis-1+BoxLineThickness; n_Y++) {
                ColorObject = new Color(
                    (int)0,
                    (int)0,
                    (int)0,
                    255
                );
                image.setRGB(n_X + MarginLeft-1 -BoxLineThickness+1, n_Y +MarginTop, ColorObject.getRGB());
                image.setRGB(n_X + MarginLeft +NumberOfPixelAlongXAxis, n_Y+MarginTop, ColorObject.getRGB());
            }
        }
        
        
        
        
        
        for (k=0;k<=NumberOfSites-1;k++) {
            image_rho[k] = new BufferedImage(
                NumberOfPixelAlongXAxis + MarginLeft + MarginRight,
                NumberOfPixelAlongYAxis + MarginTop  + MarginBottom,
                type
            );

            for(n_X = 0; n_X<=NumberOfPixelAlongXAxis-1; n_X++) {
                for(n_Y = 0; n_Y<=NumberOfPixelAlongYAxis-1; n_Y++) {
                    ColorObject = new Color(
                        (int)Math.floor(Math.random()*(double)256),
                        (int)Math.floor(Math.random()*(double)256),
                        (int)Math.floor(Math.random()*(double)256),
                        255
                    );
                    image_rho[k].setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                }
            }
            
            for(n_X = -BoxLineThickness; n_X<=NumberOfPixelAlongXAxis-1+BoxLineThickness; n_X++) {
                for(n_Y = 0; n_Y<=BoxLineThickness-1; n_Y++) {
                    ColorObject = new Color(
                        (int)0,
                        (int)0,
                        (int)0,
                        255
                    );
                    image_rho[k].setRGB(n_X + MarginLeft, n_Y +MarginTop-1 -BoxLineThickness+1, ColorObject.getRGB());
                    image_rho[k].setRGB(n_X + MarginLeft, n_Y +NumberOfPixelAlongYAxis+1 +MarginTop-1, ColorObject.getRGB());
                }
            }
            
            for(n_X = 0; n_X<=BoxLineThickness-1; n_X++) {
                for(n_Y = -BoxLineThickness; n_Y<=NumberOfPixelAlongYAxis-1+BoxLineThickness; n_Y++) {
                    ColorObject = new Color(
                        (int)0,
                        (int)0,
                        (int)0,
                        255
                    );
                    image_rho[k].setRGB(n_X + MarginLeft-1 -BoxLineThickness+1, n_Y +MarginTop, ColorObject.getRGB());
                    image_rho[k].setRGB(n_X + MarginLeft +NumberOfPixelAlongXAxis, n_Y+MarginTop, ColorObject.getRGB());
                }
            }
        }
        
        
        
        
        
        for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
            for (k=0;k<=NumberOfSites-1;k++) {
                image_absC[n][k] = new BufferedImage(
                    NumberOfPixelAlongXAxis + MarginLeft + MarginRight,
                    NumberOfPixelAlongYAxis + MarginTop  + MarginBottom,
                    type
                );

                for(n_X = 0; n_X<=NumberOfPixelAlongXAxis-1; n_X++) {
                    for(n_Y = 0; n_Y<=NumberOfPixelAlongYAxis-1; n_Y++) {
                        ColorObject = new Color(
                            (int)Math.floor(Math.random()*(double)256),
                            (int)Math.floor(Math.random()*(double)256),
                            (int)Math.floor(Math.random()*(double)256),
                            255
                        );
                        image_absC[n][k].setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                    }
                }
                
                for(n_X = -BoxLineThickness; n_X<=NumberOfPixelAlongXAxis-1+BoxLineThickness; n_X++) {
                    for(n_Y = 0; n_Y<=BoxLineThickness-1; n_Y++) {
                        ColorObject = new Color(
                            (int)0,
                            (int)0,
                            (int)0,
                            255
                        );
                        image_absC[n][k].setRGB(n_X + MarginLeft, n_Y +MarginTop-1 -BoxLineThickness+1, ColorObject.getRGB());
                        image_absC[n][k].setRGB(n_X + MarginLeft, n_Y +NumberOfPixelAlongYAxis+1 +MarginTop-1, ColorObject.getRGB());
                    }
                }
                
                for(n_X = 0; n_X<=BoxLineThickness-1; n_X++) {
                    for(n_Y = -BoxLineThickness; n_Y<=NumberOfPixelAlongYAxis-1+BoxLineThickness; n_Y++) {
                        ColorObject = new Color(
                            (int)0,
                            (int)0,
                            (int)0,
                            255
                        );
                        image_absC[n][k].setRGB(n_X + MarginLeft-1 -BoxLineThickness+1, n_Y +MarginTop, ColorObject.getRGB());
                        image_absC[n][k].setRGB(n_X + MarginLeft +NumberOfPixelAlongXAxis, n_Y+MarginTop, ColorObject.getRGB());
                    }
                }
            }
        }
        
        
        
        
        
        TimeBefore  = System.nanoTime();
        
        
        for ( n_Y=0 ; n_Y<=NumberOfPixelAlongYAxis-1; n_Y++ ) {
            mu_for = mu_max - n_Y * mu_step;

            
            
            for ( n_X=0 ; n_X<=NumberOfPixelAlongXAxis-1; n_X++ ) {
                    
                // out
                if ( (60L*60L*1000L*1000L*1000L) <= (System.nanoTime() - LastImageSaveTime) ) {
                    
                    LastImageSaveTime = System.nanoTime();
                    
                    System.out.print( "\u001B[36m\u001B[1mUpdating ./images/output/current/out.png ..." );
                    
                    try {
                        
                        ImageIO.write(image, "png",new File("./images/output/current/out "  + ".png"));
                        
                    } catch (IOException e) {
                            e.printStackTrace();
                    }
                    System.out.println( " finished\u001B[0m" );
                    
                    
                    System.out.print(
                          "\u001B[36m\u001B[1mUpdating "
                        + "./images/output/out "
                        + date.toString()
                        + ".png"
                        + " ..."
                    );
                    
                    
                    try {
                        
                        ImageIO.write(image, "png",new File("./images/output/out " + date.toString() + ".png"));
                        
                    } catch (IOException e) {
                            e.printStackTrace();
                    }
                    System.out.println( " finished\u001B[0m" );
                    
                    
                }
                
                // rho[k]
                if ( (1L*60L*1000L*1000L*1000L) <= (System.nanoTime() - LastImageSaveTimeRhos)) {
                    
                    LastImageSaveTimeRhos = System.nanoTime();
                    
                    for (k=0;k<=NumberOfSites-1;k++) {
                        System.out.print(
                              "\u001B[36m\u001B[1m"
                            + "Updating "
                            + "./images/output/current/rho" + k + ".png"
                            + " ..."
                        );
                        
                        try {
                            
                            ImageIO.write(image_rho[k], "png",new File("./images/output/current/rho" + k + ".png"));
                            
                        } catch (IOException e) {
                                e.printStackTrace();
                        }
                        System.out.println( " finished\u001B[0m" );
                    }
                    
                    for (k=0;k<=NumberOfSites-1;k++) {
                        System.out.print(
                              "\u001B[36m\u001B[1m"
                            + "Updating "
                            + "./images/output/rho" + k + " " + date.toString() + ".png"
                            + " ..."
                        );
                        try {
                            
                            ImageIO.write(image_rho[k], "png",new File("./images/output/rho" + k + " " + date.toString() + ".png"));
                            
                        } catch (IOException e) {
                                e.printStackTrace();
                        }
                        System.out.println( " finished\u001B[0m" );
                    }
                }
                
                
                // absC[k]
                if ( (20L*60L*1000L*1000L*1000L) <= (System.nanoTime() - LastImageSaveTimeAbsCs) ) {
                    
                    LastImageSaveTimeAbsCs = System.nanoTime();
                    
                    for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
                        for (k=0;k<=NumberOfSites-1;k++) {
                            System.out.print(
                                "\u001B[36m\u001B[1m"
                                + "Updating "
                                + "./images/output/current/absC_" + n + "_" + k + ".png"
                                + " ..."
                            );
                            try {
                                
                                ImageIO.write(image_absC[n][k], "png",new File("./images/output/current/absC_" + n + "_" + k + ".png"));
                                
                            } catch (IOException e) {
                                    e.printStackTrace();
                            }
                            System.out.println( " finished\u001B[0m" );
                        }
                    }
                    
                    for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
                        for (k=0;k<=NumberOfSites-1;k++) {
                            System.out.print(
                                "\u001B[36m\u001B[1m"
                                + "Updating "
                                + "./images/output/absC_" + n + "_" + k + " "  + date.toString() + ".png"
                                + " ..."
                            );
                            try {
                                
                                ImageIO.write(image_absC[n][k], "png",new File("./images/output/absC_" + n + "_" + k + " "  + date.toString() + ".png"));
                                
                            } catch (IOException e) {
                                    e.printStackTrace();
                            }
                            System.out.println( " finished\u001B[0m" );
                        }
                    }
                }
                    
                if ( (1000L*1000L*1000L) <= (System.nanoTime() - LastPrintoutTime) ) {
                    
                    LastPrintoutTime = System.nanoTime();
                    
                    TimeSoFar = System.nanoTime()-TimeBefore;
                    
                    PercentageSoFar=
                    (
                        (
                            (double)n_Y
                            + (double)n_X/(double)NumberOfPixelAlongXAxis
                        )
                        / (double)NumberOfPixelAlongYAxis
                        * 100.0
                    );
                    
                    
                    EstimatedTotalTime = TimeSoFar / (PercentageSoFar/100.0);
                    
                    EstimatedTimeLeft = EstimatedTotalTime * (1.0-PercentageSoFar / 100.0);
                    
                    System.out.print(
                        PercentageSoFar
                        + "%                          \r"
                    );
                    
                    System.out.println("");
                    
                    System.out.println(
                        "estimated time total : "
                        +(
                            EstimatedTotalTime
                            /1000.0/1000.0/1000.0/60.0
                        )
                        + " min"
                    );
                    
                    System.out.println(
                        "estimated time left  : "
                        +(
                            EstimatedTimeLeft
                            /1000.0/1000.0/1000.0/60.0
                        )
                        + " min"
                    );
                
                }
                
                J_for = J_min + n_X * J_step;
                
                EnergyFunction.SetHoppingStrength(J_for);
                EnergyFunction.SetChemicalPotential(mu_for);
                
                System.out.println("--------------------------------------------------");
                System.out.println("\u001B[35m\u001B[1mJ  = " + J_for + " | mu = " + mu_for + "\u001B[0m" );
                
                
                //StartValuesx = new double[NumberOfParameters];
                
                //StepSizesx = new double[NumberOfParameters];
                //for (n=0;n<=StepSizesx.length-1;n++) {
                //    StepSizesx[n] = 0.001;
                //}
                
                LowestEnergy=100000000.0;
                Energies = new double[NumberOfMinimizationsPerPixelInPhaseDiagram];
                
                all_goals_are_fulfilled_in_at_least_one_for_current_pixel = false;
                NumberOfSuccessesForCurrentPixel=0;
                
                CounterForTotalNumberOfMinimizationsPerPixelInPhaseDiagram=0;
                
                
                for (IterationNumber=1;IterationNumber<=NumberOfMinimizationsPerPixelInPhaseDiagram;IterationNumber++) {
                    
                    StepSizesx = new double[NumberOfParameters];
                    
                    for (n=0;n<=StepSizesx.length-1;n++) {
                        StepSizesx[n] = 0.001;
                        StepSizeModificationx[n] = StepSizeModificationDefaultValue;
                    }
                    
                    
                    for (n=0;n<=parameters.length-1;n++) {
                        parameters[n] = Math.random()*2.0-1.0;
                    }
                    
                        
                    parameters = EnergyFunction.NormalizeStateVector(parameters);
                    
                    BestSoFarInnerLoopX = parameters.clone();
                    BestSoFarInnerLoopEnergy = EnergyFunction.Value(BestSoFarInnerLoopX);
                    
                    k=0;
                    for (;;) {
                        
                        k++;
                        if (
                            Math.ceil(
                                Math.exp(
                                    a0 + a1 * (double)IterationNumber
                                )
                            )
                            < k
                        ) {
                            break;
                        };
                        
                        CounterForTotalNumberOfMinimizationsPerPixelInPhaseDiagram++;
                        
                        gradient=GradientQuadraticApproximation(EnergyFunction, parameters, StepSizesx);
                        
                        //System.out.println("-----------------------------" );
                        //System.out.println(" " + k );
                        
                        for (n=0;n<=StepSizesx.length-1;n++) {
                            
                            if ( StepSizesx[n] != 0.0) {
                                if (
                                ( Math.signum(StepSizesx[n]) != Math.signum((-1.0)*gradient[n]) )
                                ) {
                                    /* gradient changes signum */
                                    
                                    //System.out.println("Gradient changed signum for n = " + n );
                                    
                                    //System.out.print("\u001B[31m\u001B[1mStepSizeModificationx[" + n + "] changes from " + StepSizeModificationx[n] + " to " );
                                    
                                    if (
                                        StepSizeModificationMinValue <= StepSizeModificationx[n]
                                    ) {
                                        StepSizeModificationx[n]/=StepSizeModificationDecreaseValue;
                                    }
                                    
                                    if (
                                        StepSizeModificationx[n] < StepSizeModificationMinValue
                                    ) {
                                        StepSizeModificationx[n]=StepSizeModificationMinValue;
                                    }
                                    
                                    //System.out.println( StepSizeModificationx[n] + "\u001B[0m" );
                                    
                                    //System.out.println( "StepSizeModificationx[" + n + "] \u001B[32m\u001B[1mdecreased\u001B[0m" );
                                }
                                else {
                                    /* gradient doesn't change signum */
                                    
                                    //System.out.print("\u001B[32m\u001B[1mStepSizeModificationx[" + n + "] changes from " + StepSizeModificationx[n] + " to " );
                                    
                                    if (
                                        StepSizeModificationx[n] <= StepSizeModificationMaxValue
                                    ) {
                                        StepSizeModificationx[n]*=StepSizeModificationIncreaseValue;
                                    }
                                    
                                    if (
                                        StepSizeModificationMaxValue < StepSizeModificationx[n]
                                    ) {
                                        StepSizeModificationx[n]=StepSizeModificationMaxValue;
                                    }
                                    
                                    //System.out.println( StepSizeModificationx[n] + "\u001B[0m" );
                                    
                                    //System.out.println( "StepSizeModificationx[" + n + "] \u001B[31m\u001B[1mincreased\u001B[0m" );
                                }
                            
                            }
                            
                            //System.out.println( "StepSizeModificationx[" + n + "]" + StepSizeModificationx[n] );
                        }
                        
                        for (n=0;n<=StepSizesx.length-1;n++) {
                            StepSizesx[n] = (-1.0)*gradient[n]*StepSizeModificationx[n];//(0.1+Math.random()*(2.0-0.1))
                            //System.out.println("StepSizesx[" + n + "] = " + StepSizesx[n] );
                            //System.out.println("gradient[" + n + "] = " + gradient[n] );
                        }
                        
                        //System.out.println("");
                        
                        for (n=0;n<=StepSizesx.length-1;n++) {
                            parameters[n] += StepSizesx[n];
                            //parameters[n] = BestSoFarInnerLoopX[n] + StepSizesx[n];
                            //parameters[n] = LastX[n] + StepSizesx[n];
                            //System.out.println("parameters[" + n + "] = " + parameters[n] );
                        }
                        
                        
                        
                        //System.out.println("-----------------------------" );
                        //System.out.println(" " + k );
                        
                        //System.out.print( " " + k );
                        
                        //parameters = EnergyFunction.NextDatapointForMinimization(parameters);

                        parameters = EnergyFunction.NormalizeStateVector(parameters);
                        /*
                        for (n=0;n<=parameters.length-1;n++) {
                            System.out.println("parameters[" + n + "] = " + parameters[n] );
                        }*/
                        
                        //PrintPhysicalValues(parameters);
                        
                        //rho = EnergyFunction.CalculateParticleDensities(parameters);
                        
                        //for (n=0;n<=rho.length-1;n++) {
                        //    System.out.println("rho[" + n + "] = " + rho[n] );
                        //}
                        
                        //System.out.println("" );
                        
                        CurrentX = parameters.clone();
                        CurrentEnergy = EnergyFunction.Value(parameters);
                        
                        
                        //System.out.println( "k = " + k );
                        //System.out.println( "BestSoFarInnerLoopEnergy = " + BestSoFarInnerLoopEnergy );
                        //System.out.println( "CurrentEnergy = " + CurrentEnergy );
                        //System.out.println( "LastEnergy = " + LastEnergy );
                        
                        all_goals_are_fulfilled=true;
                        all_goals_are_fulfilled_involving_absolute_errors=true;
                        all_goals_are_fulfilled_involving_relative_errors=true;
                        
                        /*
                        all_goals_are_fulfilled_involving_absolute_errors &= (
                            Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy )
                            <=
                            GoalYDifferenceToLastStep
                        );
                        
                        for (n=0;n<=GoalXDifferenceToLastStep.length-1;n++) {
                            all_goals_are_fulfilled_involving_absolute_errors &= (
                                Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] )
                                <=
                                GoalXDifferenceToLastStep[n]
                            );
                        }
                        
                        all_goals_are_fulfilled_involving_relative_errors &= (
                            Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy )
                            <=
                            GoalYRelativeError * Math.abs(CurrentEnergy)
                        );
                        
                        for (n=0;n<=GoalXDifferenceToLastStep.length-1;n++) {
                            all_goals_are_fulfilled_involving_relative_errors &= (
                                Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] )
                                <=
                                GoalXRelativeError[n] * Math.abs(CurrentX[n])
                            );
                        }
                        
                        all_goals_are_fulfilled &= (
                            all_goals_are_fulfilled_involving_absolute_errors
                            |
                            all_goals_are_fulfilled_involving_relative_errors
                        );*/
                        
                        
                        /*
                        all_goals_are_fulfilled &= (
                            Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy )
                            <=
                            GoalYDifferenceToLastStep
                            +
                            GoalYRelativeError * Math.abs(CurrentEnergy)
                        );
                        
                        for (n=0;n<=GoalXDifferenceToLastStep.length-1;n++) {
                            all_goals_are_fulfilled &= (
                                Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] )
                                <=
                                GoalXDifferenceToLastStep[n]
                                +
                                GoalXRelativeError[n] * Math.abs(CurrentX[n])
                            );
                        }*/
                        
                        //System.out.println( "all_goals_are_fulfilled 1 = " + all_goals_are_fulfilled );
                        //System.out.println( "Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy ) = " + Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy ) );
                        
                        all_goals_are_fulfilled &= (
                            Math.abs( LastEnergy - CurrentEnergy )
                            <=
                            GoalYDifferenceToLastStep
                        );
                        //System.out.println( "all_goals_are_fulfilled 2 = " + all_goals_are_fulfilled );
                        
                        for (n=0;n<=GoalXDifferenceToLastStep.length-1;n++) {
                            //System.out.println( "Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] ) = " + Math.abs( LastX[n] - CurrentX[n] ) );
                            all_goals_are_fulfilled &= (
                                Math.abs( LastX[n] - CurrentX[n] )
                                <=
                                GoalXDifferenceToLastStep[n]
                            );
                        }
                        //System.out.println( "all_goals_are_fulfilled 3 = " + all_goals_are_fulfilled );
                        
                        /*
                        all_goals_are_fulfilled &= (
                            Math.abs( BestSoFarInnerLoopEnergy - CurrentEnergy )
                            <=
                            GoalYRelativeError * Math.abs(CurrentEnergy)
                        );
                        
                        for (n=0;n<=GoalXDifferenceToLastStep.length-1;n++) {
                            all_goals_are_fulfilled &= (
                                Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] )
                                <=
                                GoalXRelativeError[n] * Math.abs(CurrentX[n])
                            );
                            
                            //System.out.println( "BestSoFarInnerLoopX[n] = " + BestSoFarInnerLoopX[n] );
                            //System.out.println( "CurrentX[n] = " + CurrentX[n] );
                            //System.out.println( "Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] ) = " + Math.abs( BestSoFarInnerLoopX[n] - CurrentX[n] ) );
                            //System.out.println( "GoalXRelativeError[n] * Math.abs(CurrentX[n]) = " + GoalXRelativeError[n] * Math.abs(CurrentX[n]) );
                        }
                        */
                        
                        //System.out.println( "all_goals_are_fulfilled 3 = " + all_goals_are_fulfilled );
                        
                        
                        if (all_goals_are_fulfilled) {
                            all_goals_are_fulfilled_in_at_least_one_for_current_pixel = true;
                            NumberOfSuccessesForCurrentPixel++;
                            break;
                        }
                        
                        if ( CurrentEnergy < BestSoFarInnerLoopEnergy ) {
                            BestSoFarInnerLoopX = parameters.clone();
                            BestSoFarInnerLoopEnergy = EnergyFunction.Value(BestSoFarInnerLoopX);
                        }
                        /*
                        if ( k == 1 ) {
                            BestSoFarInnerLoopX = parameters.clone();
                            BestSoFarInnerLoopEnergy = EnergyFunction.Value(BestSoFarInnerLoopX);
                        }
                        else {
                            if ( CurrentEnergy < BestSoFarInnerLoopEnergy ) {
                                BestSoFarInnerLoopX = parameters.clone();
                                BestSoFarInnerLoopEnergy = EnergyFunction.Value(BestSoFarInnerLoopX);
                            }
                        }*/
                        
                        LastEnergy = CurrentEnergy;
                            
                        for (n=0;n<=LastX.length-1;n++) {
                            LastX[n] = CurrentX[n];
                        }
                        
                    }
                    
                    
                    
                    
                    System.out.print( " " + IterationNumber + ":" );
                    if (all_goals_are_fulfilled) {
                        System.out.print( "\u001B[32m\u001B[1m" + k + "\u001B[0m |" );
                    } else {
                        System.out.print( "\u001B[31m\u001B[1m" + k + "\u001B[0m |" );
                    }
                    
                    
                    
                    if ( IterationNumber == 1 ) {
                        BestSoFar = parameters.clone();
                    }
                    else {
                        if ( EnergyFunction.Value(BestSoFarInnerLoopX) < EnergyFunction.Value(BestSoFar) ) {
                            BestSoFar = BestSoFarInnerLoopX.clone();
                        }
                    }
                    
                    if ( NumberOfSuccessfulMinimizationsPerPixelInPhaseDiagramAfterWhichToStop <= NumberOfSuccessesForCurrentPixel ) {
                        break;
                    }
                
                }
                
                if (all_goals_are_fulfilled_in_at_least_one_for_current_pixel) {
                    System.out.print( "\nall_goals_are_fulfilled_in_at_least_one_for_current_pixel = \u001B[32m\u001B[1m " + all_goals_are_fulfilled_in_at_least_one_for_current_pixel + "\u001B[0m" );
                } else {
                    System.out.print( "\nall_goals_are_fulfilled_in_at_least_one_for_current_pixel = \u001B[31m\u001B[1m " + all_goals_are_fulfilled_in_at_least_one_for_current_pixel + "\u001B[0m" );
                }
                
                System.out.print( "\nCounterForTotalNumberOfMinimizationsPerPixelInPhaseDiagram = \u001B[36m\u001B[1m " + CounterForTotalNumberOfMinimizationsPerPixelInPhaseDiagram + "\u001B[0m" );
                
                
                System.out.println( "" );
                
                rho = EnergyFunction.CalculateParticleDensities(BestSoFar);
                
                Arrays.sort(rho);
                
                temp=EnergyFunction.Value(BestSoFar);
                
                
                OutputFileStreamWithFormatedData_Current.writeInt( n_X );
                OutputFileStreamWithFormatedData_Current.writeInt( n_Y );
                OutputFileStreamWithFormatedData_Current.writeDouble( J_for );
                OutputFileStreamWithFormatedData_Current.writeDouble( mu_for );
                OutputFileStreamWithFormatedData_Current.writeDouble( temp );
                
                OutputFileStreamWithFormatedData.writeInt( n_X );
                OutputFileStreamWithFormatedData.writeInt( n_Y );
                OutputFileStreamWithFormatedData.writeDouble( J_for );
                OutputFileStreamWithFormatedData.writeDouble( mu_for );
                OutputFileStreamWithFormatedData.writeDouble( temp );
                
                
                if ( ( ToPlotYMin <= temp) & ( temp <= ToPlotYMax ) ) {
                    red   = (int)(( temp - ToPlotYMin )/(ToPlotYMax-ToPlotYMin)*255);
                    green = red;
                    blue  = red;
                    alpha = 255;
                }
                else {
                    red   = (int)Math.floor(Math.random()*(double)256);
                    green = (int)Math.floor(Math.random()*(double)256);
                    blue  = (int)Math.floor(Math.random()*(double)256);
                    alpha = 255;
                }
                
                if ( !all_goals_are_fulfilled_in_at_least_one_for_current_pixel ) {
                    red   = 255;
                    green = 0;
                    blue  = 0;
                    alpha = 255;
                }
                
                ColorObject = new Color( red, green, blue, alpha );
                image.setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                
                /*
                InputFileStreamWithFormatedData_Current = new DataInputStream( new FileInputStream( "./output/dataset " + date.toString() + ".dat" ) );
                
                System.out.println( "" );
                
                IntReadFromFile = InputFileStreamWithFormatedData_Current.readInt();
                System.out.println( "NumberOfPixelAlongXAxis = " + IntReadFromFile );
                IntReadFromFile = InputFileStreamWithFormatedData_Current.readInt();
                System.out.println( "NumberOfPixelAlongYAxis = " + IntReadFromFile );
                System.out.println( "" );
                    
                while (true) {
                    try {
                        IntReadFromFile = InputFileStreamWithFormatedData_Current.readInt();
                        System.out.print( "( n_X = " + IntReadFromFile );
                        IntReadFromFile = InputFileStreamWithFormatedData_Current.readInt();
                        System.out.print( ", n_Y = " + IntReadFromFile );
                        DoubleReadFromFile = InputFileStreamWithFormatedData_Current.readDouble();
                        System.out.print( ", J = " + DoubleReadFromFile );
                        DoubleReadFromFile = InputFileStreamWithFormatedData_Current.readDouble();
                        System.out.print( ", mu = " + DoubleReadFromFile );
                        DoubleReadFromFile = InputFileStreamWithFormatedData_Current.readDouble();
                        System.out.print( ", E = " + DoubleReadFromFile );
                        System.out.print( " )\n" );
                    }
                    catch ( java.io.EOFException e ) {
                        InputFileStreamWithFormatedData_Current.close();
                        break;
                    }
                }
                System.out.println( "" );
                */
                
                
                
                for (k=0;k<=NumberOfSites-1;k++) {
                    temp = rho[k] ;
                    
                    if ( ( MinNumberOfParticles <= temp) & ( temp <= MaxNumberOfParticles ) ) {
                        hue        = (float)Mod(temp, 1.0);
                        saturation = 1.0f;
                        value      = 1.0f;
                    }
                    else {
                        hue        = 0.0f;
                        saturation = 0.0f;
                        value      = 0.0f;
                    }
                    
                    if ( !all_goals_are_fulfilled_in_at_least_one_for_current_pixel ) {
                        hue        = 0.0f;
                        saturation = 0.0f;
                        value      = 0.0f;
                    }
                    
                    //ColorObject = new Color( red, green, blue, alpha );
                    ColorObject = Color.getHSBColor(hue, saturation, value);
                    image_rho[k].setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                
                }
                
                for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
                    for (k=0;k<=NumberOfSites-1;k++) {
                        absC[k][n] = Math.sqrt(
                              BestSoFar[2*NumberOfCoordinatesPerSite*k+(2*n)]
                            * BestSoFar[2*NumberOfCoordinatesPerSite*k+(2*n)]
                            + BestSoFar[2*NumberOfCoordinatesPerSite*k+(2*n+1)]
                            * BestSoFar[2*NumberOfCoordinatesPerSite*k+(2*n+1)]
                        );
                    }
                }
                
                
                java.util.Arrays.sort(
                    absC,
                    new java.util.Comparator<double[]>() {
                        public int compare(double[] a, double[] b) {
                            return Double.compare(a[0], b[0]);
                        }
                    }
                );
                
                for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
                    for (k=0;k<=NumberOfSites-1;k++) {
                        temp = Mod(absC[k][n], 1.0);
                        
                        temp/= 2.0;/* using half of the color spectrum */
                        
                        if ( !all_goals_are_fulfilled_in_at_least_one_for_current_pixel ) {
                            ColorObject = Color.getHSBColor(0.0f, 0.0f, 0.0f);
                            image_absC[n][k].setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                        }
                        else {
                            ColorObject = Color.getHSBColor((float)temp, 1.0f, 1.0f);
                            image_absC[n][k].setRGB(n_X + MarginLeft, n_Y + MarginTop, ColorObject.getRGB());
                        }
                    
                    }
                }
                
                
            }
        
        }
        
        
        
        
        // out
        try {
            
            ImageIO.write(image, "png",new File("./images/output/current/out "  + ".png"));
            
        } catch (IOException e) {
                e.printStackTrace();
        }
        
        try {
            
            ImageIO.write(image, "png",new File("./images/output/out " + date.toString() + ".png"));
            
        } catch (IOException e) {
                e.printStackTrace();
        }
        
        // rho[k]
        for (k=0;k<=NumberOfSites-1;k++) {
            try {
                
                ImageIO.write(image_rho[k], "png",new File("./images/output/current/rho" + k + ".png"));
                
            } catch (IOException e) {
                    e.printStackTrace();
            }
        }
        
        for (k=0;k<=NumberOfSites-1;k++) {
            try {
                
                ImageIO.write(image_rho[k], "png",new File("./images/output/rho" + k + " " + date.toString() + ".png"));
                
            } catch (IOException e) {
                    e.printStackTrace();
            }
        }
        
        // absC[k]
        for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
            for (k=0;k<=NumberOfSites-1;k++) {
                try {
                    
                    ImageIO.write(image_absC[n][k], "png",new File("./images/output/current/absC_" + n + "_" + k + ".png"));
                    
                } catch (IOException e) {
                        e.printStackTrace();
                }
            }
        }
        
        for (n=0;n<=NumberOfCoordinatesPerSite-1;n++) {
            for (k=0;k<=NumberOfSites-1;k++) {
                try {
                    
                    ImageIO.write(image_absC[n][k], "png",new File("./images/output/absC_" + n + "_" + k + " "  + date.toString() + ".png"));
                    
                } catch (IOException e) {
                        e.printStackTrace();
                }
            }
        }
        
                
        
        TimeSoFar = System.nanoTime()-TimeBefore;
        
        PercentageSoFar=100.0;
        
        
        EstimatedTotalTime = TimeSoFar / (PercentageSoFar/100.0);
        
        EstimatedTimeLeft = EstimatedTotalTime * (1.0-PercentageSoFar / 100.0);
        
        System.out.print(
            PercentageSoFar
            + "%                          \r"
        );
        
        System.out.println("");
        
        System.out.println(
            "estimated time total : "
            +(
                EstimatedTotalTime
                /1000.0/1000.0/1000.0/60.0
            )
            + " min"
        );
        
        System.out.println(
            "estimated time left  : "
            +(
                EstimatedTimeLeft
                /1000.0/1000.0/1000.0/60.0
            )
            + " min"
        );
                
        
        System.out.println("");
        
        System.out.println( TimeSoFar + "ns" );
        System.out.println( TimeSoFar/1000.0 + "µs" );
        System.out.println( TimeSoFar/1000.0/1000.0 + "ms" );
        System.out.println( TimeSoFar/1000.0/1000.0/1000.0 + "s" );
        System.out.println( TimeSoFar/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( TimeSoFar/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( TimeSoFar/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
                        }

    }
    
    private static double func(double x) {
        return Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x)*Sin(x);
        
    }
    
    private static double[] FindLocalMinimumRoughly (double xStartValue, double StepSize, long MaxNumberOfIterations) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        long k;
        
        double x1=xStartValue;
        double x2=x1+StepSize;
        double xnew;
        double y1;
        double y2;
        double ynew;
        
        double[] toreturn = new double[3];
        
        TimeBefore = System.nanoTime();

        x1=xStartValue;
        x2=x1+StepSize;
        y1 = func(x1);
        y2 = func(x2);
        
        
        System.out.println("Iterations: " + MaxNumberOfIterations);
        for (n=1; n<=MaxNumberOfIterations; n++) {
            
            if (y1<y2) {
                xnew=x1-StepSize;
                ynew=func(xnew);
                
                if (ynew<=y1) { //going further to smaller values
                    x2=x1; x1=xnew;
                    y2=y1; y1=ynew;
                }
                else if (y1<ynew) { //local minimum found
                    System.out.println("local minimum found between " + xnew + " and " + x2 + " in step " + n );
                    toreturn[0]=xnew;
                    toreturn[1]=x1;
                    toreturn[2]=x2;
                    break;
                };
            }
            else {
                xnew=x2+StepSize;
                ynew=func(xnew);
                
                if (ynew<=y2) { //going further to bigger values
                    x1=x2; x2=xnew;
                    y1=y2; y2=ynew;
                }
                else if (y2<ynew) { //local minimum found
                    System.out.println("local minimum found between " + x1 + " and " + xnew + " in step " + n);
                    toreturn[0]=x1;
                    toreturn[1]=x2;
                    toreturn[2]=xnew;
                    break;
                };
            };
            
        };
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        return toreturn;
    }
    
    private static double[] FindLocalMinimumRoughlyWithDebuggingOutput (double xStartValue, double StepSize, long MaxNumberOfIterations) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        long k;
        
        double x1=xStartValue;
        double x2=x1+StepSize;
        double xnew;
        double y1;
        double y2;
        double ynew;
        
        double[] toreturn = new double[3];
        
        TimeBefore = System.nanoTime();

        x1=xStartValue;
        x2=x1+StepSize;
        y1 = func(x1);
        y2 = func(x2);
        
        
        System.out.println("Iterations: " + MaxNumberOfIterations);
        for (n=1; n<=MaxNumberOfIterations; n++) {
            
            if (y1<y2) {
                xnew=x1-StepSize;
                ynew=func(xnew);
                
                System.out.println(n);
                System.out.println( "(xnew, x1, x2) = (" + Round(xnew,0.0001) + ", " + Round(x1,0.0001) + ", " + Round(x2,0.0001) + ") " );
                System.out.println( "(ynew, y1, y2) = (" + Round(ynew,0.0001) + ", " + Round(y1,0.0001) + ", " + Round(y2,0.0001) + ") " );
                
                if (ynew<=y1) { //going further to smaller values
                    System.out.println("Going further to smaller values" );
                    x2=x1; x1=xnew;
                    y2=y1; y1=ynew;
                    System.out.println( "(xnew, x1, x2) = (" + Round(xnew,0.0001) + ", " + Round(x1,0.0001) + ", " + Round(x2,0.0001) + ") " );
                    System.out.println( "(ynew, y1, y2) = (" + Round(ynew,0.0001) + ", " + Round(y1,0.0001) + ", " + Round(y2,0.0001) + ") " );
                }
                else if (y1<ynew) { //local minimum found
                    System.out.println("local minimum found between " + xnew + " and " + x2 + " in step " + n );
                    toreturn[0]=xnew;
                    toreturn[1]=x1;
                    toreturn[2]=x2;
                    break;
                };
            }
            else {
                xnew=x2+StepSize;
                ynew=func(xnew);
                
                System.out.println(n);
                System.out.println( "(x1, x2, xnew) = (" + Round(x1,0.0001) + ", " + Round(x2,0.0001) + ", " + Round(xnew,0.0001) + ") " );
                System.out.println( "(y1, y2, ynew) = (" + Round(y1,0.0001) + ", " + Round(y2,0.0001) + ", " + Round(ynew,0.0001) + ") " );
                
                if (ynew<=y2) { //going further to bigger values
                    System.out.println("Going further to bigger values" );
                    x1=x2; x2=xnew;
                    y1=y2; y2=ynew;
                    System.out.println( "(x1, x2, xnew) = (" + Round(x1,0.0001) + ", " + Round(x2,0.0001) + ", " + Round(xnew,0.0001) + ") " );
                    System.out.println( "(y1, y2, ynew) = (" + Round(y1,0.0001) + ", " + Round(y2,0.0001) + ", " + Round(ynew,0.0001) + ") " );
                }
                else if (y2<ynew) { //local minimum found
                    System.out.println("local minimum found between " + x1 + " and " + xnew + " in step " + n);
                    toreturn[0]=x1;
                    toreturn[1]=x2;
                    toreturn[2]=xnew;
                    break;
                };
            };
            
            
            
            System.out.println("--------------------------------------------------");
        };
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        return toreturn;
    }

    private static double[][] SearchForLocalMinimumInIntervall(double[] SearchIntervall, int DataPointsPerIteration, long MaxNumberOfIterations, double GoalAbsoluteErrorX, double GoalRelativeErrorX, double GoalAbsoluteErrorY, double GoalRelativeErrorY, long TimeLimitInNanoSeconds) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        int k;
        
        double x1=SearchIntervall[0];
        double x2=SearchIntervall[1];
        double lastx1=SearchIntervall[0];
        double lastx2=SearchIntervall[1];
        double StepSize;
        double[][] ListOfDataPoints = new double[DataPointsPerIteration][2];
        double[][] Minimum = new double[2][2];
        
        TimeBefore = System.nanoTime();
        
        if (DataPointsPerIteration<4) {
            System.out.println(
                ">>> Number of datapoints per iteration must be at least 4.\n" +
                "    Instead found:      DataPointsPerIteration = " + DataPointsPerIteration + "\n" +
                "    Automatically set:  DataPointsPerIteration = 4\n"
            );
            DataPointsPerIteration=4;
        }
        
        for (n=1; n<=MaxNumberOfIterations; n++) {
            StepSize=
                (x2-x1)
                /
                (DataPointsPerIteration-1)
            ;
            for (k=1; k<=DataPointsPerIteration; k++) {
                ListOfDataPoints[k-1][0]=x1+(k-1)*StepSize;
                ListOfDataPoints[k-1][1]=func(ListOfDataPoints[k-1][0]);
                
                //System.out.println( ListOfDataPoints[k-1][0] + ", " + ListOfDataPoints[k-1][1] );
            }
            
            //System.out.println("------------------");
            
            Minimum=SearchForLocalMinimumInSortedListOfDataPoints(ListOfDataPoints);
            
            //System.out.println( Minimum[0][0] + ", " + Minimum[0][1]);
            //System.out.println( Minimum[1][0] + ", " + Minimum[1][1]);
            
            System.out.println( "AbsoluteErrorX = " + (Math.abs(Minimum[0][0]-Minimum[1][0])/2.0));
            System.out.println( "RelativeErrorX * Abs(x) = " + (GoalAbsoluteErrorY*Math.abs( (Minimum[0][0]+Minimum[1][0])/2.0) ));
            
            System.out.println( "AbsoluteErrorY = " + (Math.abs(Minimum[0][1]-Minimum[1][1])/2.0));
            System.out.println( "RelativeErrorY * Abs(y) = " + (GoalAbsoluteErrorY*Math.abs( (Minimum[0][1]+Minimum[1][1])/2.0) ));
            
            
            x1=Minimum[0][0];
            x2=Minimum[1][0];
            
            if (x1==lastx1 & x2==lastx2) {
                if (printing_out_debug_informations_is_activated) {
                    System.out.println("Stopped, because x-intervall doesn't change any more.");
                }
                break;
            }
            
            if (
                (
                    (
                        Math.abs(
                            (Minimum[0][0]-Minimum[1][0])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorX
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[0][0]-Minimum[1][0])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorX
                        *
                        Math.abs(
                            (Minimum[0][0]+Minimum[1][0])/2.0
                        ) 
                    )
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[0][1]-Minimum[1][1])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorY
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[0][1]-Minimum[1][1])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorY
                        *
                        Math.abs(
                            (Minimum[0][1]+Minimum[1][1])/2.0
                        )
                    )
                )
            ) {
                System.out.println("Stopped, because GoalAbsoluteErrorX, GoalRelativeErrorX, GoalAbsoluteErrorY, GoalRelativeErrorY were fulfilled.");
                System.out.println(" = " + Math.abs(
                            (Minimum[0][0]-Minimum[1][0])/2.0
                        ));
                System.out.println(" = " + Math.abs(
                            (Minimum[0][0]-Minimum[1][0])/2.0
                        )/Math.abs(
                            (Minimum[0][0]+Minimum[1][0])/2.0
                        ) );
                System.out.println(" = " + Math.abs(
                            (Minimum[0][1]-Minimum[1][1])/2.0
                        ));
                System.out.println(" = " + Math.abs(
                            (Minimum[0][1]-Minimum[1][1])/2.0
                        )/Math.abs(
                            (Minimum[0][1]+Minimum[1][1])/2.0
                        ) );
                break;
            }
            
            lastx1=x1;
            lastx2=x2;
            
            //System.out.println("------------------");
        }
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        
        return Minimum;
    }
    
    private static double[][] SearchForLocalMinimumInSortedListOfDataPoints(double[][] DataPoints) {
        double[][] toreturn = new double[2][2];
        boolean stepped_down_at_least_once = false;
        int index_of_intervall_min=0;
        int index_of_intervall_max=DataPoints.length-1;
        
        int n;
        
        for (n=1;n<=DataPoints.length-1;n++) {
            if ( DataPoints[n-1][1] < DataPoints[n][1] ) {
                index_of_intervall_max=n;
                //toreturn[1][0]=DataPoints[n][0];
                //toreturn[1][1]=DataPoints[n][1];
                if (stepped_down_at_least_once) {
                    break;
                };
            }
            else if ( DataPoints[n-1][1] > DataPoints[n][1] ) {
                stepped_down_at_least_once = true;
                index_of_intervall_min=n-1;
                //toreturn[0][0]=DataPoints[n-1][0];
                //toreturn[0][1]=DataPoints[n-1][1];
                index_of_intervall_max=n;
                //toreturn[1][0]=DataPoints[n][0];
                //toreturn[1][1]=DataPoints[n][1];
            }
        }
        
        //intervall_min
        toreturn[0][0]=DataPoints[index_of_intervall_min][0];
        toreturn[0][1]=DataPoints[index_of_intervall_max-1][1];
        //intervall_max
        toreturn[1][0]=DataPoints[index_of_intervall_max][0];
        toreturn[1][1]=Math.min(DataPoints[index_of_intervall_min][1],DataPoints[index_of_intervall_max][1]);
        
        return toreturn;
    }

    
    private static double func2(double xA, double xB) {
        return Sin(xA)*Sin(xB);
        
    }
    
    private static double[][] FindLocalMinimumRoughly2 (double StartValuex1, double StartValuex2, double StepSizex1, double StepSizex2, long MaxNumberOfIterations) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        /*
        double x1;
        double x2;
        double  y;
        
        double x1_m1_00;
        double x2_m1_00;
        double  y_m1_00;
        
        double x1_p1_00;
        double x2_p1_00;
        double  y_p1_00;
        
        double x1_00_m1;
        double x2_00_m1;
        double  y_00_m1;
        
        double x1_00_p1;
        double x2_00_p1;
        double  y_00_p1;
        */
        double[][] ListOfDataPoints = new double[5][3];
        
        double[] NewDataPoint = new double[3];
        
        double[][] toreturn = new double[3][2];
        
        TimeBefore = System.nanoTime();

        ListOfDataPoints[0][0]=StartValuex1;
        ListOfDataPoints[0][1]=StartValuex2;
        ListOfDataPoints[0][2]=func2(ListOfDataPoints[0][0],ListOfDataPoints[0][1]);
        
        ListOfDataPoints[1][0]=ListOfDataPoints[0][0]-StepSizex1;
        ListOfDataPoints[1][1]=ListOfDataPoints[0][1];
        ListOfDataPoints[1][2]=func2(ListOfDataPoints[1][0],ListOfDataPoints[1][1]);
        
        ListOfDataPoints[2][0]=ListOfDataPoints[0][0]+StepSizex1;
        ListOfDataPoints[2][1]=ListOfDataPoints[0][1];
        ListOfDataPoints[2][2]=func2(ListOfDataPoints[2][0],ListOfDataPoints[2][1]);
        
        ListOfDataPoints[3][0]=ListOfDataPoints[0][0];
        ListOfDataPoints[3][1]=ListOfDataPoints[0][1]-StepSizex2;
        ListOfDataPoints[3][2]=func2(ListOfDataPoints[3][0],ListOfDataPoints[3][1]);
        
        ListOfDataPoints[4][0]=ListOfDataPoints[0][0];
        ListOfDataPoints[4][1]=ListOfDataPoints[0][1]+StepSizex2;
        ListOfDataPoints[4][2]=func2(ListOfDataPoints[4][0],ListOfDataPoints[4][1]);
        
        
        System.out.println("Iterations: " + MaxNumberOfIterations);
        
        for (n=1; n<=MaxNumberOfIterations; n++) {
            /*System.out.println("-----------------------------------------------");
            System.out.println("n=" + n);
            
            System.out.println("[" + ListOfDataPoints[0][0] + ", " + ListOfDataPoints[0][1] + ", " + ListOfDataPoints[0][2] + "] ");
            System.out.println("");
            
            System.out.println("[" + ListOfDataPoints[1][0] + ", " + ListOfDataPoints[1][1] + ", " + ListOfDataPoints[1][2] + "] ");
            System.out.println("[" + ListOfDataPoints[2][0] + ", " + ListOfDataPoints[2][1] + ", " + ListOfDataPoints[2][2] + "] ");
            System.out.println("");
            
            System.out.println("[" + ListOfDataPoints[3][0] + ", " + ListOfDataPoints[3][1] + ", " + ListOfDataPoints[3][2] + "] ");
            System.out.println("[" + ListOfDataPoints[4][0] + ", " + ListOfDataPoints[4][1] + ", " + ListOfDataPoints[4][2] + "] ");
            System.out.println("");
            */
            
            NewDataPoint=FindDatapointWithLowestValueInListOfDataPoints(ListOfDataPoints);
            
            
            /*System.out.println("NewDataPoint:");
            System.out.println("[" + NewDataPoint[0] + ", " + NewDataPoint[1] + ", " + NewDataPoint[2] + "] ");
            System.out.println("");
            */
            
            if (
            (NewDataPoint[0] == ListOfDataPoints[0][0]) &
            (NewDataPoint[1] == ListOfDataPoints[0][1]) &
            (NewDataPoint[2] == ListOfDataPoints[0][2])
            ) {
            
                break;
            }
            
            
            ListOfDataPoints[0][0]=NewDataPoint[0];
            ListOfDataPoints[0][1]=NewDataPoint[1];
            ListOfDataPoints[0][2]=NewDataPoint[2];
            
            ListOfDataPoints[1][0]=ListOfDataPoints[0][0]-StepSizex1;
            ListOfDataPoints[1][1]=ListOfDataPoints[0][1];
            ListOfDataPoints[1][2]=func2(ListOfDataPoints[1][0],ListOfDataPoints[1][1]);
            
            ListOfDataPoints[2][0]=ListOfDataPoints[0][0]+StepSizex1;
            ListOfDataPoints[2][1]=ListOfDataPoints[0][1];
            ListOfDataPoints[2][2]=func2(ListOfDataPoints[2][0],ListOfDataPoints[2][1]);
            
            ListOfDataPoints[3][0]=ListOfDataPoints[0][0];
            ListOfDataPoints[3][1]=ListOfDataPoints[0][1]-StepSizex2;
            ListOfDataPoints[3][2]=func2(ListOfDataPoints[3][0],ListOfDataPoints[3][1]);
            
            ListOfDataPoints[4][0]=ListOfDataPoints[0][0];
            ListOfDataPoints[4][1]=ListOfDataPoints[0][1]+StepSizex2;
            ListOfDataPoints[4][2]=func2(ListOfDataPoints[4][0],ListOfDataPoints[4][1]);
            
        };
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println("finished after " + n + " loops");
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        double[] ListOfYValues={
            ListOfDataPoints[1][2],
            ListOfDataPoints[2][2],
            ListOfDataPoints[3][2],
            ListOfDataPoints[4][2]
        };
        
        // y intervall
        toreturn[0][0]=NewDataPoint[2]; /* lower boundary */
        toreturn[0][1]=Max(ListOfYValues); /* upper boundary */
        
        // x1 intervall
        toreturn[1][0]=ListOfDataPoints[1][0]; /* lower boundary */
        toreturn[1][1]=ListOfDataPoints[2][0]; /* upper boundary */
        
        // x2 intervall
        toreturn[2][0]=ListOfDataPoints[3][1]; /* lower boundary */
        toreturn[2][1]=ListOfDataPoints[4][1]; /* upper boundary */
        
        return toreturn;
    }
    
    
    private static double[] FindDatapointWithLowestValueInListOfDataPoints(double[][] ListOfDataPoints) {
        
        double[] CurrentlyLargest = new double[3];
        
        int n;
        
        CurrentlyLargest[0] = ListOfDataPoints[0][0];
        CurrentlyLargest[1] = ListOfDataPoints[0][1];
        CurrentlyLargest[2] = ListOfDataPoints[0][2];
        
        for (n=1; n<=ListOfDataPoints.length-1; n++) {
            if ( ListOfDataPoints[n][2] < CurrentlyLargest[2] ) {
                CurrentlyLargest[0] = ListOfDataPoints[n][0];
                CurrentlyLargest[1] = ListOfDataPoints[n][1];
                CurrentlyLargest[2] = ListOfDataPoints[n][2];
                //System.out.println("\n[" + CurrentlyLargest[0] + ", " + CurrentlyLargest[1] + ", " + CurrentlyLargest[2] + "] ");
            }
            else {
                //System.out.println("\n[" + CurrentlyLargest[0] + ", " + CurrentlyLargest[1] + ", " + CurrentlyLargest[2] + "] ");
            }
        }
        
        return CurrentlyLargest;
        
    }
    
    private static double[][] SearchForLocalMinimumInIntervall2(double[][] SearchIntervall, int DataPointsPerIterationPerDimension, long MaxNumberOfIterations, double GoalAbsoluteErrorX1, double GoalRelativeErrorX1, double GoalAbsoluteErrorX2, double GoalRelativeErrorX2, double GoalAbsoluteErrorY, double GoalRelativeErrorY, long TimeLimitInNanoSeconds) {
        
        long TimeBefore;
        long TimeAfter;
        
        long IterationNumber;
        int DatapointNumberInX1Direction;
        int DatapointNumberInX2Direction;
        int Direction;
        
        double x11=SearchIntervall[1][0];
        double x12=SearchIntervall[1][1];
        double lastx11=SearchIntervall[1][0];
        double lastx12=SearchIntervall[1][1];
        double x21=SearchIntervall[2][0];
        double x22=SearchIntervall[2][1];
        double lastx21=SearchIntervall[2][0];
        double lastx22=SearchIntervall[2][1];
        double StepSizeInX1Direction;
        double StepSizeInX2Direction;
        double[][][] ListOfDataPoints = new double[2][DataPointsPerIterationPerDimension][3];
        double[][][] Minimum = new double[3][2][2];
        
        TimeBefore = System.nanoTime();
        
        if (DataPointsPerIterationPerDimension<4) {
            System.out.println(
                ">>> Number of datapoints per iteration must be at least 4.\n" +
                "    Instead found:      DataPointsPerIterationPerDimension = " + DataPointsPerIterationPerDimension + "\n" +
                "    Automatically set:  DataPointsPerIterationPerDimension = 4\n"
            );
            DataPointsPerIterationPerDimension=4;
        }
        
        for (IterationNumber=1; IterationNumber<=MaxNumberOfIterations; IterationNumber++) {
            StepSizeInX1Direction=
                (x12-x11)
                /
                (DataPointsPerIterationPerDimension-1)
            ;
            StepSizeInX2Direction=
                (x22-x21)
                /
                (DataPointsPerIterationPerDimension-1)
            ;
            
            /* Calculating datapoints in x1 direction */
            
            double x2ValueForDatapointsInX1Direction = (x22+x21)/2;
            
            double x1;
            
            for (DatapointNumberInX1Direction=1; DatapointNumberInX1Direction<=DataPointsPerIterationPerDimension; DatapointNumberInX1Direction++) {
                
                /* current x1 value */
                x1=x11 + (DatapointNumberInX1Direction-1) * StepSizeInX1Direction;
                //System.out.println("x1 = " + x1);
                
                /* y values for datapoint list in x1 direction */
                ListOfDataPoints[1-1][DatapointNumberInX1Direction-1][0]=func2(
                        x1,
                        x2ValueForDatapointsInX1Direction
                );
                
                /* x1 values for datapoint list in x1 direction */
                ListOfDataPoints[1-1][DatapointNumberInX1Direction-1][1]=x1;
                
                /* x2 value for datapoint list in x1 direction */
                ListOfDataPoints[1-1][DatapointNumberInX1Direction-1][2]=x2ValueForDatapointsInX1Direction;
                
            }
            
            
            /* Calculating datapoints in x2 direction */
            
            double x1ValueForDatapointsInX2Direction = (x12+x11)/2;
            
            double x2;
            
            for (DatapointNumberInX2Direction=1; DatapointNumberInX2Direction<=DataPointsPerIterationPerDimension; DatapointNumberInX2Direction++) {
                
                /* current x2 value */
                x2=x21 + (DatapointNumberInX2Direction-1) * StepSizeInX2Direction;
                //System.out.println("x2 = " + x2);
                
                /* y values for datapoint list in x2 direction */
                ListOfDataPoints[2-1][DatapointNumberInX2Direction-1][0]=func2(
                        x1ValueForDatapointsInX2Direction,
                        x2
                );
                
                /* x1 values for datapoint list in x2 direction */
                ListOfDataPoints[2-1][DatapointNumberInX2Direction-1][1]=x1ValueForDatapointsInX2Direction;
                
                /* x2 value for datapoint list in x2 direction */
                ListOfDataPoints[2-1][DatapointNumberInX2Direction-1][2]=x2;
                
            }
            
            
            
            for (Direction=1;Direction<=2;Direction++) {
                Minimum[Direction]=SearchForLocalMinimumInSortedListOfDataPoints2(ListOfDataPoints,Direction);
            }
            
            double CurrentY=Minimum[1][0][0];
            
            /* find min y */
            double CurrentlySmallestY = Minimum[1][0][0];
            for (Direction=1;Direction<=2;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][0][0];
                }
                else {
                    CurrentY=Minimum[Direction][1][0];
                }
                
                if ( CurrentY < CurrentlySmallestY ) {
                    CurrentlySmallestY = CurrentY;
                }
            }
            Minimum[0][0][0]=CurrentY;
            
            
            /* find max y */
            double CurrentlyLargestY = Minimum[1][0][0];
            for (Direction=1;Direction<=2;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][1][0];
                }
                else {
                    CurrentY=Minimum[Direction][0][0];
                }
                
                if ( CurrentY < CurrentlyLargestY ) {
                    CurrentlyLargestY = CurrentY;
                }
            }
            Minimum[0][1][0]=CurrentY;
            
            
            
            x11=Minimum[1][0][1];
            x12=Minimum[1][1][1];
            
            x21=Minimum[2][0][1];
            x22=Minimum[2][1][1];
            
            
            if (x11==lastx11 & x12==lastx12 & x21==lastx21 & x22==lastx22) {
                System.out.println("Stopped, because x_n-intervalls don't change any more.");
                break;
            }
            
            if (
                (
                    (
                        Math.abs(
                            (Minimum[1][0][1]-Minimum[1][1][1])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorX1
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[1][0][1]-Minimum[1][1][1])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorX1
                        *
                        Math.abs(
                            (Minimum[1][0][1]+Minimum[1][1][1])/2.0
                        ) 
                    )
                ) &
                
                (
                    (
                        Math.abs(
                            (Minimum[1][0][0]-Minimum[1][1][0])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorY
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[1][0][0]-Minimum[1][1][0])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorY
                        *
                        Math.abs(
                            (Minimum[1][0][0]+Minimum[1][1][0])/2.0
                        )
                    )
                ) &
                
                (
                    (
                        Math.abs(
                            (Minimum[2][0][1]-Minimum[2][1][1])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorX1
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[2][0][1]-Minimum[2][1][1])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorX1
                        *
                        Math.abs(
                            (Minimum[2][0][1]+Minimum[2][1][1])/2.0
                        ) 
                    )
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[2][0][0]-Minimum[2][1][0])/2.0
                        )
                    )
                    <
                    GoalAbsoluteErrorY
                ) &
                (
                    (
                        Math.abs(
                            (Minimum[2][0][0]-Minimum[2][1][0])/2.0
                        )
                    )
                    <
                    (
                        GoalRelativeErrorY
                        *
                        Math.abs(
                            (Minimum[2][0][0]+Minimum[2][1][0])/2.0
                        )
                    )
                )

            ) {
                
                break;
            }
            
            lastx11=Minimum[1][0][1];
            lastx12=Minimum[1][1][1];
            
            lastx21=Minimum[2][0][1];
            lastx22=Minimum[2][1][1];
        }
        
        
        double[][] toreturn = new double[3][2];
        
        toreturn[0][0] = Minimum[0][0][0];
        toreturn[0][1] = Minimum[0][1][0];
            
        for (Direction=1;Direction<=2;Direction++) {
            toreturn[Direction][0] = Minimum[Direction][0][1];
            toreturn[Direction][1] = Minimum[Direction][1][1];
        }
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        
        return toreturn;
    }

    
    private static double[][] SearchForLocalMinimumInSortedListOfDataPoints2(double[][][] DataPoints, int Direction) {
        double[][] toreturn = new double[2][2];
        boolean stepped_down_at_least_once = false;
        int index_of_intervall_min=0;
        int index_of_intervall_max=DataPoints[Direction-1].length-1;
        
        int n;
        
        for (n=1;n<=DataPoints[Direction-1].length-1;n++) {
            if ( DataPoints[Direction-1][n-1][0] < DataPoints[Direction-1][n][0] ) {
                index_of_intervall_max=n;
                if (stepped_down_at_least_once) {
                    break;
                };
            }
            else if ( DataPoints[Direction-1][n-1][0] > DataPoints[Direction-1][n][0] ) {
                stepped_down_at_least_once = true;
                index_of_intervall_min=n-1;
                index_of_intervall_max=n;
            }
        }
        
        //intervall_min
        toreturn[0][0]=DataPoints[Direction-1][index_of_intervall_max-1][0];
        toreturn[0][1]=DataPoints[Direction-1][index_of_intervall_min][Direction];
        //intervall_max
        toreturn[1][0]=Math.min(DataPoints[Direction-1][index_of_intervall_min][0],DataPoints[Direction-1][index_of_intervall_max][0]);
        toreturn[1][1]=DataPoints[Direction-1][index_of_intervall_max][Direction];
        
        return toreturn;
    }
    
    private static double funcN(double[] x) {
        
        double toreturn=0.0;
        int n;
        
        for (n=0;n<=x.length-1;n++) {
            toreturn += Sin(x[n]);
            //System.out.println("n = " + n);
        }
        
        return toreturn;
        
    }
    
    
    private static double[][] FindLocalMinimumRoughlyN (double[] StartValuesx, double[] StepSizesx, long MaxNumberOfIterations) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        int k;
        int Direction;
        
        double[][] ListOfDataPoints = new double[2*StartValuesx.length + 1][StartValuesx.length + 1];
        
        double[] NewDataPoint = new double[StartValuesx.length + 1];
        
        double[][] toreturn = new double[StartValuesx.length + 1][2];
        
        TimeBefore = System.nanoTime();
        
        double[] parameters = new double[StartValuesx.length];
        
        for (Direction=1;Direction<=StartValuesx.length;Direction++) {
            ListOfDataPoints[0][Direction]=StartValuesx[Direction-1];
            parameters[Direction-1]=ListOfDataPoints[0][Direction];
        }
        ListOfDataPoints[0][0]=funcN(parameters);
        /*
        System.out.println("ListOfDataPoints[0][0] = " + ListOfDataPoints[0][0]);
        System.out.println("ListOfDataPoints[0][1] = " + ListOfDataPoints[0][1]);
        System.out.println("ListOfDataPoints[0][2] = " + ListOfDataPoints[0][2]);
        System.out.println("ListOfDataPoints[0][3] = " + ListOfDataPoints[0][3]);
        System.out.println("");
        */
        
        for (k=1;k<=StartValuesx.length;k++) {
        
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                if (Direction==k) {
                    ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction]-StepSizesx[k-1];
                }
                else {
                    ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction];
                }
                parameters[Direction-1]=ListOfDataPoints[2*k-1][Direction];
            }
            ListOfDataPoints[2*k-1][0]=funcN(parameters);
            /*
            System.out.println("StepSizesx[" + (k-1) + "] = " + StepSizesx[k-1]);
            System.out.println("ListOfDataPoints[" + (2*k-1) + "][0] = " + ListOfDataPoints[2*k-1][0]);
            System.out.println("ListOfDataPoints[" + (2*k-1) + "][1] = " + ListOfDataPoints[2*k-1][1]);
            System.out.println("ListOfDataPoints[" + (2*k-1) + "][2] = " + ListOfDataPoints[2*k-1][2]);
            System.out.println("ListOfDataPoints[" + (2*k-1) + "][3] = " + ListOfDataPoints[2*k-1][3]);
            System.out.println("");
            */
            
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                if (Direction==k) {
                    ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction]+StepSizesx[k-1];
                }
                else {
                    ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction];
                }
                parameters[Direction-1]=ListOfDataPoints[1][Direction];
            }
            ListOfDataPoints[2*k][0]=funcN(parameters);
            /*
            System.out.println("StepSizesx[" + (k-1) + "] = " + StepSizesx[k-1]);
            System.out.println("ListOfDataPoints[" + (2*k) + "][0] = " + ListOfDataPoints[2*k][0]);
            System.out.println("ListOfDataPoints[" + (2*k) + "][1] = " + ListOfDataPoints[2*k][1]);
            System.out.println("ListOfDataPoints[" + (2*k) + "][2] = " + ListOfDataPoints[2*k][2]);
            System.out.println("ListOfDataPoints[" + (2*k) + "][3] = " + ListOfDataPoints[2*k][3]);
            System.out.println("");
            */
        }
        
        System.out.println("Iterations: " + MaxNumberOfIterations);
        
        for (n=1; n<=MaxNumberOfIterations-1; n++) {
            /*
            System.out.println("-----------------------------------------------");
            System.out.println("n=" + n);
            
            System.out.println("[" + ListOfDataPoints[0][0] + ", " + ListOfDataPoints[0][1] + ", " + ListOfDataPoints[0][2] + "] ");
            System.out.println("");
            
            System.out.println("[" + ListOfDataPoints[1][0] + ", " + ListOfDataPoints[1][1] + ", " + ListOfDataPoints[1][2] + "] ");
            System.out.println("[" + ListOfDataPoints[2][0] + ", " + ListOfDataPoints[2][1] + ", " + ListOfDataPoints[2][2] + "] ");
            System.out.println("");
            
            System.out.println("[" + ListOfDataPoints[3][0] + ", " + ListOfDataPoints[3][1] + ", " + ListOfDataPoints[3][2] + "] ");
            System.out.println("[" + ListOfDataPoints[4][0] + ", " + ListOfDataPoints[4][1] + ", " + ListOfDataPoints[4][2] + "] ");
            System.out.println("");
            */
            
            NewDataPoint=FindDatapointWithLowestValueInListOfDataPointsN(ListOfDataPoints);
            
            
            /*System.out.println("NewDataPoint:");
            System.out.println("[" + NewDataPoint[0] + ", " + NewDataPoint[1] + ", " + NewDataPoint[2] + "] ");
            System.out.println("");
            */
            
            boolean minimum_is_found=true;
            
            for (Direction=0;Direction<=StartValuesx.length;Direction++) {
                minimum_is_found &= ( NewDataPoint[Direction] == ListOfDataPoints[0][Direction] );
            }
            
            if (minimum_is_found) {break;}
            
            for (Direction=0;Direction<=StartValuesx.length;Direction++) {
                ListOfDataPoints[0][Direction]=NewDataPoint[Direction];
            }
            /*
            System.out.println("ListOfDataPoints[0][0] = " + ListOfDataPoints[0][0]);
            System.out.println("ListOfDataPoints[0][1] = " + ListOfDataPoints[0][1]);
            System.out.println("ListOfDataPoints[0][2] = " + ListOfDataPoints[0][2]);
            System.out.println("ListOfDataPoints[0][3] = " + ListOfDataPoints[0][3]);
            System.out.println("");
            */
            
            for (k=1;k<=StartValuesx.length;k++) {
            
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    if (Direction==k) {
                        ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction]-StepSizesx[k-1];
                    }
                    else {
                        ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction];
                    }
                    parameters[Direction-1]=ListOfDataPoints[2*k-1][Direction];
                }
                ListOfDataPoints[2*k-1][0]=funcN(parameters);
                /*
                System.out.println("StepSizesx[" + (k-1) + "] = " + StepSizesx[k-1]);
                System.out.println("ListOfDataPoints[" + (2*k-1) + "][0] = " + ListOfDataPoints[2*k-1][0]);
                System.out.println("ListOfDataPoints[" + (2*k-1) + "][1] = " + ListOfDataPoints[2*k-1][1]);
                System.out.println("ListOfDataPoints[" + (2*k-1) + "][2] = " + ListOfDataPoints[2*k-1][2]);
                System.out.println("ListOfDataPoints[" + (2*k-1) + "][3] = " + ListOfDataPoints[2*k-1][3]);
                System.out.println("");
                */
                
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    if (Direction==k) {
                        ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction]+StepSizesx[k-1];
                    }
                    else {
                        ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction];
                    }
                    parameters[Direction-1]=ListOfDataPoints[1][Direction];
                }
                ListOfDataPoints[2*k][0]=funcN(parameters);
                /*
                System.out.println("StepSizesx[" + (k-1) + "] = " + StepSizesx[k-1]);
                System.out.println("ListOfDataPoints[" + (2*k) + "][0] = " + ListOfDataPoints[2*k][0]);
                System.out.println("ListOfDataPoints[" + (2*k) + "][1] = " + ListOfDataPoints[2*k][1]);
                System.out.println("ListOfDataPoints[" + (2*k) + "][2] = " + ListOfDataPoints[2*k][2]);
                System.out.println("ListOfDataPoints[" + (2*k) + "][3] = " + ListOfDataPoints[2*k][3]);
                System.out.println("");
                */
            }
            /*
            ListOfDataPoints[0][0]=NewDataPoint[0];
            ListOfDataPoints[0][1]=NewDataPoint[1];
            ListOfDataPoints[0][2]=NewDataPoint[2];
            
            ListOfDataPoints[1][1]=ListOfDataPoints[0][1]-StepSizex1;
            ListOfDataPoints[1][2]=ListOfDataPoints[0][2];
            ListOfDataPoints[1][0]=func2(ListOfDataPoints[1][1],ListOfDataPoints[1][2]);
            
            ListOfDataPoints[2][1]=ListOfDataPoints[0][1]+StepSizex1;
            ListOfDataPoints[2][2]=ListOfDataPoints[0][2];
            ListOfDataPoints[2][0]=func2(ListOfDataPoints[2][2],ListOfDataPoints[2][2]);
            
            ListOfDataPoints[3][1]=ListOfDataPoints[0][1];
            ListOfDataPoints[3][2]=ListOfDataPoints[0][2]-StepSizex2;
            ListOfDataPoints[3][0]=func2(ListOfDataPoints[3][1],ListOfDataPoints[3][2]);
            
            ListOfDataPoints[4][1]=ListOfDataPoints[0][1];
            ListOfDataPoints[4][2]=ListOfDataPoints[0][2]+StepSizex2;
            ListOfDataPoints[4][0]=func2(ListOfDataPoints[4][1],ListOfDataPoints[4][2]);
            */
        };
        
        TimeAfter = System.nanoTime();
        
        System.out.println("");
        
        System.out.println("finished after " + n + " loops");
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        double[] ListOfYValues = new double[StartValuesx.length];
        
        for (k=1;k<=StartValuesx.length-1;k++) {
            ListOfYValues[k]=ListOfDataPoints[k][0];
        }
        
        // y intervall
        toreturn[0][0]=NewDataPoint[0]; /* lower boundary */
        toreturn[0][1]=Max(ListOfYValues); /* upper boundary */
        
        // x_n intervalls
        for (k=1;k<=StartValuesx.length;k++) {
            toreturn[k][0]=ListOfDataPoints[2*k-1][k]; /* lower boundary */
            toreturn[k][1]=ListOfDataPoints[2*k][k]; /* upper boundary */
        }
        
        return toreturn;
    }
    
    
    private static double[] FindDatapointWithLowestValueInListOfDataPointsN(double[][] ListOfDataPoints) {
        
        double[] CurrentlyLargest = new double[ListOfDataPoints[0].length];
        
        int n;
        int k;
        
        for (k=0;k<=CurrentlyLargest.length-1;k++) {
            CurrentlyLargest[k] = ListOfDataPoints[0][k];
        }
        
        for (n=1; n<=ListOfDataPoints.length-1; n++) {
            if ( ListOfDataPoints[n][0] < CurrentlyLargest[0] ) {
                for (k=0;k<=CurrentlyLargest.length-1;k++) {
                    CurrentlyLargest[k] = ListOfDataPoints[n][k];
                }
            }
        }
        
        return CurrentlyLargest;
        
    }
    
    
    private static double[][] SearchForLocalMinimumInIntervallN(
    double[][] SearchIntervall,
    int DataPointsPerIterationPerDimension,
    long MaxNumberOfIterations,
    double[] GoalAbsoluteError,
    double[] GoalRelativeError,
    double GoalAbsoluteErrorY,
    double GoalRelativeErrorY,
    long TimeLimitInNanoSeconds
    ) {
        
        long TimeBefore;
        long TimeAfter;
        
        long IterationNumber;
        int DatapointNumber;
        int Direction;
        
        /*
        double x11=SearchIntervall[1][0];
        double x12=SearchIntervall[1][1];
        double lastx11=SearchIntervall[1][0];
        double lastx12=SearchIntervall[1][1];
        double x21=SearchIntervall[2][0];
        double x22=SearchIntervall[2][1];
        double lastx21=SearchIntervall[2][0];
        double lastx22=SearchIntervall[2][1];
        double StepSizeInX1Direction;
        double StepSizeInX2Direction;
        */
        
        double[] xA = new double[SearchIntervall.length];
        double[] xB = new double[SearchIntervall.length];
        double[] lastxA = new double[SearchIntervall.length];
        double[] lastxB = new double[SearchIntervall.length];
        double[] StepSize = new double[SearchIntervall.length];
        
        for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
            xA[Direction]=SearchIntervall[Direction][0];
            xB[Direction]=SearchIntervall[Direction][1];
            lastxA[Direction]=SearchIntervall[Direction][0];
            lastxB[Direction]=SearchIntervall[Direction][1];
        }
        
        double[][][] ListOfDataPoints = new double[SearchIntervall.length][DataPointsPerIterationPerDimension][SearchIntervall.length+1];
        
        double[][][] Minimum = new double[SearchIntervall.length+1][2][2];
        
        TimeBefore = System.nanoTime();
        
        if (DataPointsPerIterationPerDimension<4) {
            System.out.println(
                ">>> Number of datapoints per iteration must be at least 4.\n" +
                "    Instead found:      DataPointsPerIterationPerDimension = " + DataPointsPerIterationPerDimension + "\n" +
                "    Automatically set:  DataPointsPerIterationPerDimension = 4\n"
            );
            DataPointsPerIterationPerDimension=4;
        }
        
        for (IterationNumber=1; IterationNumber<=MaxNumberOfIterations; IterationNumber++) {
            /*StepSizeInX1Direction=
                (x12-x11)
                /
                (DataPointsPerIterationPerDimension-1)
            ;
            StepSizeInX2Direction=
                (x22-x21)
                /
                (DataPointsPerIterationPerDimension-1)
            ;*/
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                StepSize[Direction]=
                (xB[Direction]-xA[Direction])
                /
                (DataPointsPerIterationPerDimension-1);
            }
            
            /* Calculating datapoints */
            
            //double x2ValueForDatapointsInX1Direction = (x22+x21)/2;
            double[] AverageXValues = new double[SearchIntervall.length];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                AverageXValues[Direction] = (xB[Direction]+xA[Direction])/2.0;
            }
            
            double[] x = new double[SearchIntervall.length];
            double[] parameters;
            int Direction2;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                for (DatapointNumber=1; DatapointNumber<=DataPointsPerIterationPerDimension; DatapointNumber++) {
                    
                    /* current x[Direction] value */
                    x[Direction]=xA[Direction] + (DatapointNumber-1) * StepSize[Direction];
                    
                    /* y values for datapoint list in x_Direction direction */
                    parameters=AverageXValues.clone();
                    parameters[Direction-1]=x[Direction];
                    ListOfDataPoints[Direction-1][DatapointNumber-1][0]=funcN(parameters);
                    
                    /* x_n values for datapoint list in x_Direction direction */
                    for (Direction2=1;Direction2<=SearchIntervall.length-1;Direction2++) {
                        if (Direction==Direction2) {
                            ListOfDataPoints[Direction-1][DatapointNumber-1][Direction2]=x[Direction];
                        }
                        else{
                            ListOfDataPoints[Direction-1][DatapointNumber-1][Direction2]=AverageXValues[Direction2];
                        }
                    }
                    
                }
            }
            
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                Minimum[Direction]=SearchForLocalMinimumInSortedListOfDataPointsN(ListOfDataPoints,Direction);
            }
            
            double CurrentY=Minimum[1][0][0];
            
            /* find min y */
            double CurrentlySmallestY = Minimum[1][0][0];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][0][0];
                }
                else {
                    CurrentY=Minimum[Direction][1][0];
                }
                
                if ( CurrentY < CurrentlySmallestY ) {
                    CurrentlySmallestY = CurrentY;
                }
            }
            Minimum[0][0][0]=CurrentY;
            
            
            /* find max y */
            double CurrentlyLargestY = Minimum[1][0][0];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][1][0];
                }
                else {
                    CurrentY=Minimum[Direction][0][0];
                }
                
                if ( CurrentY < CurrentlyLargestY ) {
                    CurrentlyLargestY = CurrentY;
                }
            }
            Minimum[0][1][0]=CurrentY;
            
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                xA[Direction]=Minimum[Direction][0][1];
                xB[Direction]=Minimum[Direction][1][1];
            }
            
            boolean x_intervalls_do_not_change_any_more=true;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                x_intervalls_do_not_change_any_more &= (
                    xA[Direction] == lastxA[Direction] &
                    xB[Direction] == lastxB[Direction]
                );
            }
            
            if (x_intervalls_do_not_change_any_more) {
                System.out.println("Stopped, because x_n-intervalls don't change any more.");
                break;
            }
            
            boolean all_goals_are_fulfilled=true;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                all_goals_are_fulfilled &= (
                
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][1]-Minimum[Direction][1][1])/2.0
                            )
                        )
                        <
                        GoalAbsoluteError[Direction]
                    ) &
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][1]-Minimum[Direction][1][1])/2.0
                            )
                        )
                        <
                        (
                            GoalRelativeError[Direction]
                            *
                            Math.abs(
                                (Minimum[Direction][0][1]+Minimum[Direction][1][1])/2.0
                            ) 
                        )
                    ) &
                    
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][0]-Minimum[Direction][1][0])/2.0
                            )
                        )
                        <
                        GoalAbsoluteErrorY
                    ) &
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][0]-Minimum[Direction][1][0])/2.0
                            )
                        )
                        <
                        (
                            GoalRelativeErrorY
                            *
                            Math.abs(
                                (Minimum[Direction][0][0]+Minimum[Direction][1][0])/2.0
                            )
                        )
                    )
                    
                );
            }
            
            if (all_goals_are_fulfilled) {break;}
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                lastxA[Direction]=Minimum[Direction][0][1];
                lastxB[Direction]=Minimum[Direction][1][1];
            }
        }
        
        double[][] toreturn = new double[SearchIntervall.length+1][2];
        
        toreturn[0][0] = Minimum[0][0][0];
        toreturn[0][1] = Minimum[0][1][0];
            
        for (Direction=1;Direction<=SearchIntervall.length;Direction++) {
            toreturn[Direction][0] = Minimum[Direction][0][1];
            toreturn[Direction][1] = Minimum[Direction][1][1];
        }
        
        TimeAfter = System.nanoTime();

        System.out.println("loops:" + IterationNumber);
                
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        
        
        return toreturn;
    }
    
    
    private static double[][] SearchForLocalMinimumInSortedListOfDataPointsN(double[][][] DataPoints, int Direction) {
        double[][] toreturn = new double[2][2];
        boolean stepped_down_at_least_once = false;
        int index_of_intervall_min=0;
        int index_of_intervall_max=DataPoints[Direction-1].length-1;
        
        int n;
        
        for (n=1;n<=DataPoints[Direction-1].length-1;n++) {
            if ( DataPoints[Direction-1][n-1][0] < DataPoints[Direction-1][n][0] ) {
                index_of_intervall_max=n;
                if (stepped_down_at_least_once) {
                    break;
                };
            }
            else if ( DataPoints[Direction-1][n-1][0] > DataPoints[Direction-1][n][0] ) {
                stepped_down_at_least_once = true;
                index_of_intervall_min=n-1;
                index_of_intervall_max=n;
            }
        }
        
        //intervall_min
        toreturn[0][0]=DataPoints[Direction-1][index_of_intervall_max-1][0];
        toreturn[0][1]=DataPoints[Direction-1][index_of_intervall_min][Direction];
        //intervall_max
        toreturn[1][0]=Math.min(DataPoints[Direction-1][index_of_intervall_min][0],DataPoints[Direction-1][index_of_intervall_max][0]);
        toreturn[1][1]=DataPoints[Direction-1][index_of_intervall_max][Direction];
        
        return toreturn;
    }
    
    
    private static double energy001(double[] x, double[] parameter) {
        
        double toreturn=0.0;
        int n;
        
        
        double J  = parameter[0];
        double V  = parameter[1];
        double mu = parameter[2];
        
        double[] rho = new double[2];
        double[] Repsi = new double[2];
        double[] Impsi = new double[2];
        
        double[] ReC0 = new double[2];
        double[] ImC0 = new double[2];
        double[] ReC1 = new double[2];
        double[] ImC1 = new double[2];
        
        double[][] hopping = new double[2][2];
        
        double NormSquared;
        
        ReC0[0] = x[0]; ImC0[0] = x[1];
        ReC1[0] = x[2]; ImC1[0] = x[3];
        
        ReC0[1] = x[4]; ImC0[1] = x[5];
        ReC1[1] = x[6]; ImC1[1] = x[7];
        
        
        NormSquared = ReC0[0]*ReC0[0] + ImC0[0]*ImC0[0] +
                      ReC1[0]*ReC1[0] + ImC1[0]*ImC1[0] +
                      ReC0[1]*ReC0[1] + ImC0[1]*ImC0[1] +
                      ReC1[1]*ReC1[1] + ImC1[1]*ImC1[1]
        ;
        
        rho[0]=(ReC1[0]*ReC1[0]+ImC1[0]*ImC1[0]) /NormSquared;
        rho[1]=(ReC1[1]*ReC1[1]+ImC1[1]*ImC1[1]) /NormSquared;
        
        /*
        psi[0]  = (ReC0[0]-i ImC0[0])(ReC1[0]+i ImC1[0])
                = (ReC0[0] ReC1[0]-i ImC0[0] ReC1[0]+ReC0[0] i ImC1[0]-i ImC0[0] i ImC1[0])
                = (ReC0[0] ReC1[0]-i ImC0[0] ReC1[0]+ReC0[0] i ImC1[0]+ ImC0[0] ImC1[0])
        */
        
        Repsi[0] = (ReC0[0]*ReC1[0] + ImC0[0]*ImC1[0]) /NormSquared;
        Impsi[0] = (ReC0[0]*ImC1[0] - ImC0[0]*ReC1[0]) /NormSquared;
        
        Repsi[1] = (ReC0[1]*ReC1[1] + ImC0[1]*ImC1[1]) /NormSquared;
        Impsi[1] = (ReC0[1]*ImC1[1] - ImC0[1]*ReC1[1]) /NormSquared;
        
        /*
        hopping[0][1]   = psic[0]*psi[1]+psic[1]*psi[0]
                        = psic[0]*psi[1]+
                          psic[1]*psi[0]
                        = (Repsi[0]-i Impsi[0])*(Repsi[1]+i Impsi[1])+
                          (Repsi[1]-i Impsi[1])*(Repsi[0]+i Impsi[0])
                        = (Repsi[0] Repsi[1]-i Impsi[0] Repsi[1]+Repsi[0] i Impsi[1]-i Impsi[0]i Impsi[1])+
                          (Repsi[1] Repsi[0]-i Impsi[1] Repsi[0]+Repsi[1] i Impsi[0]-i Impsi[1] i Impsi[0])
                        = Repsi[0] Repsi[1]-i Impsi[0] Repsi[1]+i Repsi[0] Impsi[1]+ Impsi[0] Impsi[1]+
                          Repsi[1] Repsi[0]-i Impsi[1] Repsi[0]+i Repsi[1] Impsi[0]+Impsi[1] Impsi[0]
                        = Repsi[0] Repsi[1]+ Impsi[0] Impsi[1]+Repsi[1] Repsi[0]+Impsi[1] Impsi[0]+
                          -i Impsi[0] Repsi[1]+i Repsi[0] Impsi[1]-i Impsi[1] Repsi[0]+i Repsi[1] Impsi[0]
                        = 2 Repsi[0] Repsi[1]+2 Impsi[0] Impsi[1]
        */
        
        hopping[0][1] = 2.0*(Repsi[0]*Repsi[1]+Impsi[0]*Impsi[1]);
        
        toreturn=
            -J * (
                hopping[0][1]
            )
            +V/2.0*(
                2.0*rho[0]*rho[1]
            )
            -mu * (
                rho[0]+
                rho[1]
            )
        ;
        
        return toreturn;
        
    }
    
    
    private static double[][] FindLocalMinimumRoughly001 (double[] StartValuesx, double[] StepSizesx, long MaxNumberOfIterations, Energy EnergyFunction) {
        long TimeBefore;
        long TimeAfter;
        
        long n;
        int k;
        int Direction;
        
        double[][] ListOfDataPoints = new double[2*StartValuesx.length + 1][StartValuesx.length + 1];
        
        double[] NewDataPoint = new double[StartValuesx.length + 1];
        
        double[][] toreturn = new double[StartValuesx.length + 1][2];
        
        TimeBefore = System.nanoTime();
        
        double[] parameters = new double[StartValuesx.length];
        
        for (Direction=1;Direction<=StartValuesx.length;Direction++) {
            ListOfDataPoints[0][Direction]=StartValuesx[Direction-1];
            parameters[Direction-1]=ListOfDataPoints[0][Direction];
        }
        parameters = EnergyFunction.NormalizeStateVector(parameters);
        for (Direction=1;Direction<=StartValuesx.length;Direction++) {
            ListOfDataPoints[0][Direction]=parameters[Direction-1];
        }
        ListOfDataPoints[0][0]=EnergyFunction.Value(parameters);
        
        for (k=1;k<=StartValuesx.length;k++) {
        
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                if (Direction==k) {
                    ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction]-StepSizesx[k-1];
                }
                else {
                    ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction];
                }
                parameters[Direction-1]=ListOfDataPoints[2*k-1][Direction];
            }
            parameters = EnergyFunction.NormalizeStateVector(parameters);
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                ListOfDataPoints[2*k-1][Direction]=parameters[Direction-1];
            }
            ListOfDataPoints[2*k-1][0]=EnergyFunction.Value(parameters);
            
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                if (Direction==k) {
                    ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction]+StepSizesx[k-1];
                }
                else {
                    ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction];
                }
                parameters[Direction-1]=ListOfDataPoints[1][Direction];
            }
            parameters = EnergyFunction.NormalizeStateVector(parameters);
            for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                ListOfDataPoints[2*k][Direction]=parameters[Direction-1];
            }
            ListOfDataPoints[2*k][0]=EnergyFunction.Value(parameters);
        }
        
        //System.out.println("Iterations: " + MaxNumberOfIterations);
        
        for (n=1; n<=MaxNumberOfIterations-1; n++) {
            
            NewDataPoint=FindDatapointWithLowestValueInListOfDataPointsN(ListOfDataPoints);
            
            boolean minimum_is_found=true;
            
            for (Direction=0;Direction<=StartValuesx.length;Direction++) {
                minimum_is_found &= ( NewDataPoint[Direction] == ListOfDataPoints[0][Direction] );
            }
            
            if (minimum_is_found) {break;}
            
            for (Direction=0;Direction<=StartValuesx.length;Direction++) {
                ListOfDataPoints[0][Direction]=NewDataPoint[Direction];
            }
            
            for (k=1;k<=StartValuesx.length;k++) {
            
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    if (Direction==k) {
                        ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction]-StepSizesx[k-1];
                    }
                    else {
                        ListOfDataPoints[2*k-1][Direction]=ListOfDataPoints[0][Direction];
                    }
                    parameters[Direction-1]=ListOfDataPoints[2*k-1][Direction];
                }
                parameters = EnergyFunction.NormalizeStateVector(parameters);
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    ListOfDataPoints[2*k-1][Direction]=parameters[Direction-1];
                }
                ListOfDataPoints[2*k-1][0]=EnergyFunction.Value(parameters);
                
                
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    if (Direction==k) {
                        ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction]+StepSizesx[k-1];
                    }
                    else {
                        ListOfDataPoints[2*k][Direction]=ListOfDataPoints[0][Direction];
                    }
                    parameters[Direction-1]=ListOfDataPoints[1][Direction];
                }
                parameters = EnergyFunction.NormalizeStateVector(parameters);
                for (Direction=1;Direction<=StartValuesx.length;Direction++) {
                    ListOfDataPoints[2*k][Direction]=parameters[Direction-1];
                }
                ListOfDataPoints[2*k][0]=EnergyFunction.Value(parameters);
            }
        };
        
        if (printing_out_debug_informations_is_activated) {
            if ( MaxNumberOfIterations <= n ) {
                System.out.println( "\u001B[36m\u001B[1mFindLocalMinimumRoughly001\u001B[0m: \u001B[31m\u001B[1mnot successful\u001B[0m after " + n + " loops" );
            }
            else {
                System.out.println( "\u001B[36m\u001B[1mFindLocalMinimumRoughly001\u001B[0m: \u001B[32m\u001B[1msuccessful\u001B[0m after " + n + " loops" );
            }
        }
        
        TimeAfter = System.nanoTime();
        /*
        System.out.println("");
        
        System.out.println("finished after " + n + " loops");
        
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        */
        double[] ListOfYValues = new double[StartValuesx.length];
        
        for (k=1;k<=StartValuesx.length-1;k++) {
            ListOfYValues[k]=ListOfDataPoints[k][0];
        }
        
        // y intervall
        toreturn[0][0]=NewDataPoint[0]; /* lower boundary */
        toreturn[0][1]=Max(ListOfYValues); /* upper boundary */
        
        // x_n intervalls
        for (k=1;k<=StartValuesx.length;k++) {
            toreturn[k][0]=ListOfDataPoints[2*k-1][k]; /* lower boundary */
            toreturn[k][1]=ListOfDataPoints[2*k  ][k]; /* upper boundary */
        }
        
        return toreturn;
    }
    
    
    private static double[][] SearchForLocalMinimumInIntervall001(
    double[][] SearchIntervall,
    int DataPointsPerIterationPerDimension,
    long MaxNumberOfIterations,
    double[] GoalAbsoluteError,
    double[] GoalRelativeError,
    double GoalAbsoluteErrorY,
    double GoalRelativeErrorY,
    long TimeLimitInNanoSeconds,
    Energy EnergyFunction
    ) {
        
        long TimeBefore;
        long TimeAfter;
        
        long IterationNumber;
        int DatapointNumber;
        int Direction;
        
        /*
        double x11=SearchIntervall[1][0];
        double x12=SearchIntervall[1][1];
        double lastx11=SearchIntervall[1][0];
        double lastx12=SearchIntervall[1][1];
        double x21=SearchIntervall[2][0];
        double x22=SearchIntervall[2][1];
        double lastx21=SearchIntervall[2][0];
        double lastx22=SearchIntervall[2][1];
        double StepSizeInX1Direction;
        double StepSizeInX2Direction;
        */
        
        double[] xA = new double[SearchIntervall.length-1];
        double[] xB = new double[SearchIntervall.length-1];
        double[] lastxA = new double[SearchIntervall.length-1];
        double[] lastxB = new double[SearchIntervall.length-1];
        double[] StepSize = new double[SearchIntervall.length-1];
        
        for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
            xA[Direction-1]=SearchIntervall[Direction][0];
            xB[Direction-1]=SearchIntervall[Direction][1];
            lastxA[Direction-1]=SearchIntervall[Direction][0];
            lastxB[Direction-1]=SearchIntervall[Direction][1];
                //System.out.println("xA[" + (Direction-1) + "] = " + xA[Direction-1]);
                //System.out.println("xB[" + (Direction-1) + "] = " + xB[Direction-1]);
        }
        
        double[][][] ListOfDataPoints = new double[SearchIntervall.length][DataPointsPerIterationPerDimension][SearchIntervall.length+1];
        
        double[][][] Minimum = new double[SearchIntervall.length+1][2][2];
        
        TimeBefore = System.nanoTime();
        
        if (DataPointsPerIterationPerDimension<4) {
            System.out.println(
                ">>> Number of datapoints per iteration must be at least 4.\n" +
                "    Instead found:      DataPointsPerIterationPerDimension = " + DataPointsPerIterationPerDimension + "\n" +
                "    Automatically set:  DataPointsPerIterationPerDimension = 4\n"
            );
            DataPointsPerIterationPerDimension=4;
        }
        
        for (IterationNumber=1; IterationNumber<=MaxNumberOfIterations; IterationNumber++) {
            /*StepSizeInX1Direction=
                (x12-x11)
                /
                (DataPointsPerIterationPerDimension-1)
            ;
            StepSizeInX2Direction=
                (x22-x21)
                /
                (DataPointsPerIterationPerDimension-1)
            ;*/
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                StepSize[Direction-1]=
                (xB[Direction-1]-xA[Direction-1])
                /
                (DataPointsPerIterationPerDimension-1);
            }
            
            /* Calculating datapoints */
            
            //double x2ValueForDatapointsInX1Direction = (x22+x21)/2;
            double[] AverageXValues = new double[SearchIntervall.length-1];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                AverageXValues[Direction-1] = (xB[Direction-1]+xA[Direction-1])/2.0;
                //System.out.println("xA[" + (Direction-1) + "] = " + xA[Direction-1]);
                //System.out.println("xB[" + (Direction-1) + "] = " + xB[Direction-1]);
                //System.out.println("AverageXValues[" + (Direction-1) + "] = " + AverageXValues[Direction-1]);
            }
            //System.out.println("................");
            
            double[] x = new double[SearchIntervall.length-1];
            double[] parameters;
            int Direction2;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                for (DatapointNumber=1; DatapointNumber<=DataPointsPerIterationPerDimension; DatapointNumber++) {
                    
                    /* current x[Direction] value */
                    x[Direction-1]=xA[Direction-1] + (DatapointNumber-1) * StepSize[Direction-1];
                    
                    /* y values for datapoint list in x_Direction direction */
                    parameters=AverageXValues.clone();
                    parameters[Direction-1]=x[Direction-1];
                    parameters=EnergyFunction.NormalizeStateVector(parameters);
                    //PrintPhysicalValues(parameters);
                    ListOfDataPoints[Direction-1][DatapointNumber-1][0]=EnergyFunction.Value(parameters);
                    
                    /* x_n values for datapoint list in x_Direction direction */
                    for (Direction2=1;Direction2<=SearchIntervall.length-1;Direction2++) {
                        //if (Direction==Direction2) {
                        //    ListOfDataPoints[Direction-1][DatapointNumber-1][Direction2]=x[Direction];
                        //}
                        //else{
                        //    ListOfDataPoints[Direction-1][DatapointNumber-1][Direction2]=AverageXValues[Direction2];
                            ListOfDataPoints[Direction-1][DatapointNumber-1][Direction2]=parameters[Direction2-1];
                        //}
                    }
                    
                }
            }
            
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                Minimum[Direction]=SearchForLocalMinimumInSortedListOfDataPointsN(ListOfDataPoints,Direction);
            }
            
            double CurrentY=Minimum[1][0][0];
            
            /* find min y */
            double CurrentlySmallestY = Minimum[1][0][0];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][0][0];
                }
                else {
                    CurrentY=Minimum[Direction][1][0];
                }
                
                if ( CurrentY < CurrentlySmallestY ) {
                    CurrentlySmallestY = CurrentY;
                }
            }
            Minimum[0][0][0]=CurrentY;
            
            
            /* find max y */
            double CurrentlyLargestY = Minimum[1][0][0];
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                
                if (Minimum[Direction][0][0] < Minimum[Direction][1][0]) {
                    CurrentY=Minimum[Direction][1][0];
                }
                else {
                    CurrentY=Minimum[Direction][0][0];
                }
                
                if ( CurrentY < CurrentlyLargestY ) {
                    CurrentlyLargestY = CurrentY;
                }
            }
            Minimum[0][1][0]=CurrentY;
            
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                xA[Direction-1]=Minimum[Direction][0][1];
                xB[Direction-1]=Minimum[Direction][1][1];
            }
            
            boolean x_intervalls_do_not_change_any_more=true;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                x_intervalls_do_not_change_any_more &= (
                    xA[Direction-1] == lastxA[Direction-1] &
                    xB[Direction-1] == lastxB[Direction-1]
                );
            }
            
            if (x_intervalls_do_not_change_any_more) {
                if (printing_out_debug_informations_is_activated) {
                    System.out.println( "\u001B[32m\u001B[1m" + System.nanoTime() + "\u001B[0m Stopped, because x_n-intervalls don't change any more.");
                }
                break;
            }
            
            boolean all_goals_are_fulfilled=true;
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                all_goals_are_fulfilled &= (
                
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][1]-Minimum[Direction][1][1])/2.0
                            )
                        )
                        <
                        GoalAbsoluteError[Direction]
                    ) &
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][1]-Minimum[Direction][1][1])/2.0
                            )
                        )
                        <
                        (
                            GoalRelativeError[Direction]
                            *
                            Math.abs(
                                (Minimum[Direction][0][1]+Minimum[Direction][1][1])/2.0
                            ) 
                        )
                    ) &
                    
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][0]-Minimum[Direction][1][0])/2.0
                            )
                        )
                        <
                        GoalAbsoluteErrorY
                    ) &
                    (
                        (
                            Math.abs(
                                (Minimum[Direction][0][0]-Minimum[Direction][1][0])/2.0
                            )
                        )
                        <
                        (
                            GoalRelativeErrorY
                            *
                            Math.abs(
                                (Minimum[Direction][0][0]+Minimum[Direction][1][0])/2.0
                            )
                        )
                    )
                    
                );
            }
            
            if (all_goals_are_fulfilled) {break;}
            
            for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
                lastxA[Direction-1]=Minimum[Direction][0][1];
                lastxB[Direction-1]=Minimum[Direction][1][1];
            }
        }
        
        if (printing_out_debug_informations_is_activated) {
            if ( MaxNumberOfIterations <= IterationNumber ) {
                System.out.println( "\u001B[36m\u001B[1mSearchForLocalMinimumInIntervall001\u001B[0m: \u001B[31m\u001B[1mnot successful\u001B[0m after " + IterationNumber + " loops" );
            }
            else {
                System.out.println( "\u001B[36m\u001B[1mSearchForLocalMinimumInIntervall001\u001B[0m: \u001B[32m\u001B[1msuccessful\u001B[0m after " + IterationNumber + " loops" );
            }
        }
        
        double[][] toreturn = new double[SearchIntervall.length][2];
        
        toreturn[0][0] = Minimum[0][0][0];
        toreturn[0][1] = Minimum[0][1][0];
            
        for (Direction=1;Direction<=SearchIntervall.length-1;Direction++) {
            toreturn[Direction][0] = Minimum[Direction][0][1];
            toreturn[Direction][1] = Minimum[Direction][1][1];
        }
        
        TimeAfter = System.nanoTime();
        /*
        System.out.println("loops:" + IterationNumber);
                
        System.out.println("");
        
        System.out.println( (TimeAfter-TimeBefore) + "ns" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0 + "µs" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0 + "ms" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0 + "s" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0 + " min" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0 + " h" );
        System.out.println( (TimeAfter-TimeBefore)/1000.0/1000.0/1000.0/60.0/60.0/24.0 + " d" );
        System.out.println("--------------------------------------------------");
        */
        
        return toreturn;
    }
    
    private static double[] NormalizeStateVector_old (double[] StateVector) {
        //double NormSquared=0.0;
        double[] NormSquared = new double[2];
        int n;
        
        NormSquared[0]=0.0;
        
        for (n=0;n<=3;n++) {
            NormSquared[0] += StateVector[n]*StateVector[n];
        }
        
        NormSquared[1]=0.0;
        
        for (n=4;n<=7;n++) {
            NormSquared[1] += StateVector[n]*StateVector[n];
        }
        
        for (n=0;n<=3;n++) {
            StateVector[n] /= Math.sqrt(NormSquared[0]);
            //System.out.println("StateVector[" + n + "]" + StateVector[n]);
        }
        
        for (n=4;n<=7;n++) {
            StateVector[n] /= Math.sqrt(NormSquared[1]);
            //System.out.println("StateVector[" + n + "]" + StateVector[n]);
        }
        
        //System.out.println("oooooooooooooooooooooooooooooooooooo");
        
        return StateVector;
        
    }
    
    
    private static void PrintNormOfStateVector (double[] StateVector) {
        //double NormSquared=0.0;
        double[] NormSquared = new double[2];
        int n;
        
        NormSquared[0]=0.0;
        
        for (n=0;n<=3;n++) {
            NormSquared[0] += StateVector[n]*StateVector[n];
        }
        
        NormSquared[1]=0.0;
        
        for (n=4;n<=7;n++) {
            NormSquared[1] += StateVector[n]*StateVector[n];
        }
        
        
        System.out.println(">>> NormSquared[0] = " + NormSquared[0] );
        System.out.println(">>> NormSquared[1] = " + NormSquared[1] );
    }
    
    
    private static void PrintPhysicalValues(double[] x) {
        
        int NumberOfSites=x.length/4;
        int n;
        
        double[] rho = new double[NumberOfSites];
        double[] Repsi = new double[NumberOfSites];
        double[] Impsi = new double[NumberOfSites];
        
        double[] ReC0 = new double[NumberOfSites];
        double[] ImC0 = new double[NumberOfSites];
        double[] ReC1 = new double[NumberOfSites];
        double[] ImC1 = new double[NumberOfSites];
        
            double[] NormSquared = new double[2];
            
            ReC0[0] = x[0]; ImC0[0] = x[1];
            ReC1[0] = x[2]; ImC1[0] = x[3];
            
            ReC0[1] = x[4]; ImC0[1] = x[5];
            ReC1[1] = x[6]; ImC1[1] = x[7];
            
            
            NormSquared[0] =    ReC0[0]*ReC0[0] + ImC0[0]*ImC0[0] +
                                ReC1[0]*ReC1[0] + ImC1[0]*ImC1[0]
            ;
            
            NormSquared[1] =    ReC0[1]*ReC0[1] + ImC0[1]*ImC0[1] +
                                ReC1[1]*ReC1[1] + ImC1[1]*ImC1[1]
            ;
            
            rho[0]=(ReC1[0]*ReC1[0]+ImC1[0]*ImC1[0]) /NormSquared[0];
            rho[1]=(ReC1[1]*ReC1[1]+ImC1[1]*ImC1[1]) /NormSquared[1];
            
            /*
            psi[0]  = (ReC0[0]-i ImC0[0])(ReC1[0]+i ImC1[0])
                    = (ReC0[0] ReC1[0]-i ImC0[0] ReC1[0]+ReC0[0] i ImC1[0]-i ImC0[0] i ImC1[0])
                    = (ReC0[0] ReC1[0]-i ImC0[0] ReC1[0]+ReC0[0] i ImC1[0]+ ImC0[0] ImC1[0])
            */
            
            Repsi[0] = (ReC0[0]*ReC1[0] + ImC0[0]*ImC1[0]) /NormSquared[0];
            Impsi[0] = (ReC0[0]*ImC1[0] - ImC0[0]*ReC1[0]) /NormSquared[0];
            
            Repsi[1] = (ReC0[1]*ReC1[1] + ImC0[1]*ImC1[1]) /NormSquared[1];
            Impsi[1] = (ReC0[1]*ImC1[1] - ImC0[1]*ReC1[1]) /NormSquared[1];
            
        System.out.println( "");
        for (n=1;n<=x.length;n++) {
            System.out.println( "x[" + n + "] = " + x[n-1] );
        }
        
        System.out.println( "");
        System.out.println( "C0[0] = " + ReC0[0] + "+i " + ImC0[0] );
        System.out.println( "C1[0] = " + ReC1[0] + "+i " + ImC1[0] );
        System.out.println( "C0[1] = " + ReC0[1] + "+i " + ImC0[1] );
        System.out.println( "C1[1] = " + ReC1[1] + "+i " + ImC1[1] );
        
        System.out.println( "");
        System.out.format("rho[0] = %1.4f%n", rho[0]);
        System.out.format("rho[1] = %1.4f%n", rho[1]);
        System.out.println( "rho[0]-rho[1] = " + (rho[0]-rho[1]) );
        System.out.println( "(rho[0]-rho[1]) / (rho[0]+rho[1]) = " + ( (rho[0]-rho[1]) / (rho[0]+rho[1]) ) );
        
        System.out.println( "");
        System.out.println( "psi[0] = " + Repsi[0] + "+i " + Impsi[0] );
        System.out.println( "psi[1] = " + Repsi[1] + "+i " + Impsi[1] );
        System.out.println( "Arg(psi[0])/(2*pi) = " + Math.atan2(Repsi[0], Impsi[0])/Math.PI );
        System.out.println( "Arg(psi[1])/(2*pi) = " + Math.atan2(Repsi[1], Impsi[1])/Math.PI );
        System.out.format("(Arg(psi[0])-Arg(psi[1]))/(2*pi) = %1.4f%n", (  (Math.atan2(Repsi[0], Impsi[0])-Math.atan2(Repsi[1], Impsi[1]))/Math.PI  ));
        
        System.out.println( "");
        System.out.println( "Norm[0] = " + Math.sqrt(NormSquared[0]) );
        System.out.println( "Norm[1] = " + Math.sqrt(NormSquared[1]) );
        
        System.out.println( "");
    }
    
    private static double[] GradientLinearApproximation(Energy EnergyFunction, double[] x, double[] StepSize) {
        
        double[] gradient = new double[x.length];
        
        double[] xMiddle;
        double[] xPlus;
        double[] xMinus;
        
        //double y;
        double yMiddle;
        double yPlus;
        double yMinus;
        
        int n;
        
        //y = EnergyFunction.Value(x);
        
        
        xMiddle = x.clone();
        yMiddle = EnergyFunction.Value(xMiddle);
        
        for (n=0;n<=x.length-1;n++) {
            
            xPlus = x.clone();
            //xPlus[n] += StepSize[n];
            xPlus[n] += Math.pow(10.0, -10);//Math.pow(10.0, -10)
            yPlus = EnergyFunction.Value(xPlus);
            
            xMinus = x.clone();
            //xMinus[n] -= StepSize[n];
            xMinus[n] -= Math.pow(10.0, -10);//Math.pow(10.0, -10)
            yMinus = EnergyFunction.Value(xMinus);
            /*
            if (
                  ( xPlus[n] - xMinus[n] ) != 0
            ) {
                gradient[n] =
                    ( yPlus - yMinus )/( xPlus[n] - xMinus[n] )
                ;
            }
            else {
                gradient[n] = 0.0;
            }
            */
            
            if (
                  ( ( xMiddle[n] - xMinus[n] ) != 0.0 )
                & ( ( xPlus[n]   - xMiddle[n] ) != 0.0 )
            ) {
                gradient[n] =
                    (
                        ( yMiddle - yMinus )/( xMiddle[n] - xMinus[n] )
                    )
                    +
                    (
                        ( yPlus   - yMiddle )/( xPlus[n]  - xMiddle[n] )
                    )
                    /2.0
                ;
            }
            else {
                gradient[n] = 0.0;
            }
            
        }
        
        return gradient;
    }
    
    
    private static double[] GradientQuadraticApproximation(Energy EnergyFunction, double[] x, double[] StepSize) {
        
        double[] gradient = new double[x.length];
        
        double[] xMiddle;
        double[] xPlus;
        double[] xMinus;
        
        double xDifference12;
        double xDifference23;
        double xDifference31;
        
        
        //double xDifferenceOfSquares12;
        //double xDifferenceOfSquares23;
        //double xDifferenceOfSquares31;
        
        double yMiddle;
        double yPlus;
        double yMinus;
        
        double Numerator;
        double Denominator;
        
        int n;
        
        //y = EnergyFunction.Value(x);
        
        for (n=0;n<=x.length-1;n++) {
        
            xMiddle = x.clone();
            yMiddle = EnergyFunction.Value(x);
            
            xPlus = x.clone();
            xPlus[n] += Math.pow(10.0, -10);//Math.pow(10.0, -10)
            yPlus = EnergyFunction.Value(xPlus);
            
            xMinus = x.clone();
            xMinus[n] -= Math.pow(10.0, -10);//Math.pow(10.0, -10)
            yMinus = EnergyFunction.Value(xMinus);
            
            xDifference12 = xMinus[n]  - xMiddle[n] ;
            xDifference23 = xMiddle[n] - xPlus[n]   ;
            xDifference31 = xPlus[n]   - xMinus[n]  ;
            
            //xDifferenceOfSquares12 = xMinus[n]*xMinus[n]   - xMiddle[n]*xMiddle[n] ;
            //xDifferenceOfSquares23 = xMiddle[n]*xMiddle[n] - xPlus[n]*xPlus[n]     ;
            //xDifferenceOfSquares31 = xPlus[n]*xPlus[n]     - xMinus[n]*xMinus[n]   ;
            
            
            Denominator = (
                  xDifference12
                * xDifference23
                * xDifference31
            );
            
            Numerator = (
                  yMinus  * (-1)*xDifference23*xDifference23
                + yMiddle * ( xDifference23*xDifference23 - xDifference12*xDifference12 )
                + yPlus   * xDifference12*xDifference12
            );
            
            
            if (
                  Denominator != 0.0
            ) {
                gradient[n] = (
                    Numerator/Denominator
                );
                //System.out.println( "Denominator = " + Denominator );
                //System.out.println( "Numerator = " + Numerator );
                //System.out.println( "gradient[" + n + "] = " + gradient[n] );
                //System.out.println( "" );
            }
            else {
                gradient[n] = 0.0;
            }
            
        }
        
        return gradient;
    }
    
    
}



    class Energy
    {
        private double J;
        private double V;
        private double mu;
        
        private int Nsites;
        private int MinNumberOfParticles;
        private int MaxNumberOfParticles;
        private int NumberOfCoordinatesPerSite;
        private int NumberOfRealValuesPerSite;
        
        private int[] HoppingDistributionMatrixIndex1;
        private int[] HoppingDistributionMatrixIndex2;
        
        private double[] HoppingDistributionMatrixValue_RealPart;
        private double[] HoppingDistributionMatrixValue_ImaginaryPart;
        
        //private double toreturn=0.0;
        private int n;
        private int m;
        
        private double[] rho;
        private double[] Repsi;
        private double[] Impsi;
        private double[] onsite;
        
        private double[] Abspsi;
        /*
        private double[] ReC0;
        private double[] ImC0;
        private double[] ReC1;
        private double[] ImC1;
        */
        private double[][] ReC;
        private double[][] ImC;
        
        private double[] NormSquared;
            
        
        public Energy(int arg_MinNumberOfParticles, int arg_MaxNumberOfParticles) {
            J  = 0.0;
            V  = 0.0;
            mu = 0.0;
            
            Nsites=3;
            
            MinNumberOfParticles=arg_MinNumberOfParticles;
            MaxNumberOfParticles=arg_MaxNumberOfParticles;
            
            NumberOfCoordinatesPerSite = MaxNumberOfParticles-MinNumberOfParticles+1;
            NumberOfRealValuesPerSite = 2*NumberOfCoordinatesPerSite;
            
            rho    = new double[Nsites];
            Repsi  = new double[Nsites];
            Impsi  = new double[Nsites];
            onsite = new double[Nsites];
            
            Abspsi = new double[Nsites];
            /*
            ReC0 = new double[Nsites];
            ImC0 = new double[Nsites];
            ReC1 = new double[Nsites];
            ImC1 = new double[Nsites];
            */
            ReC = new double[NumberOfCoordinatesPerSite][Nsites];
            ImC = new double[NumberOfCoordinatesPerSite][Nsites];
            
            NormSquared = new double[Nsites];
            
            HoppingDistributionMatrixIndex1 = new int[3];
            HoppingDistributionMatrixIndex2 = new int[3];
            
            HoppingDistributionMatrixValue_RealPart      = new double[3];
            HoppingDistributionMatrixValue_ImaginaryPart = new double[3];
            
            HoppingDistributionMatrixIndex1[0]              = 1-1;
            HoppingDistributionMatrixIndex2[0]              = 2-1;
            HoppingDistributionMatrixValue_RealPart[0]      = 0.5;
            //HoppingDistributionMatrixValue_ImaginaryPart[0] = 0.0;
            
            HoppingDistributionMatrixIndex1[1]              = 2-1;
            HoppingDistributionMatrixIndex2[1]              = 3-1;
            HoppingDistributionMatrixValue_RealPart[1]      = 0.5;
            //HoppingDistributionMatrixValue_ImaginaryPart[1] = 0.0;
            
            HoppingDistributionMatrixIndex1[2]              = 3-1;
            HoppingDistributionMatrixIndex2[2]              = 1-1;
            HoppingDistributionMatrixValue_RealPart[2]      = 0.5;
            //HoppingDistributionMatrixValue_ImaginaryPart[2] = 0.0;
        }
        
        public void SetHoppingStrength(double HoppingStrength) {
            J=HoppingStrength;
        }
        
        public double HoppingStrength() {
            return J;
        }
        
        public void SetInterSiteInteractionStrength(double InterSiteInteraction) {
            V=InterSiteInteraction;
        }
        
        public double InterSiteInteractionStrength() {
            return V;
        }
        
        public void SetChemicalPotential(double ChemicalPotential) {
            mu=ChemicalPotential;
        }
        
        public double ChemicalPotential() {
            return mu;
        }
        
        public double Value(double[] Coordinates) {
            return this.energy(Coordinates);
        }
        
        
        public double[] NormalizeStateVector (double[] StateVector) {
            //double NormSquared=0.0;
            double[] NormSquared = new double[Nsites];
            int n;
            int m;
            /*
            NormSquared[0]=0.0;
            
            for (n=0;n<=3;n++) {
                NormSquared[0] += StateVector[n]*StateVector[n];
            }
            
            NormSquared[1]=0.0;
            
            for (n=4;n<=7;n++) {
                NormSquared[1] += StateVector[n]*StateVector[n];
            }
            */
            
            for (n=0;n<=NormSquared.length-1;n++) {
                NormSquared[n]=0.0;
                
                for (m=n*NumberOfRealValuesPerSite;m<=(n+1)*NumberOfRealValuesPerSite-1;m++) {
                    NormSquared[n] += StateVector[m]*StateVector[m];
                }
            }
            
            
            /*
            for (n=0;n<=3;n++) {
                StateVector[n] /= Math.sqrt(NormSquared[0]);
                //System.out.println("StateVector[" + n + "]" + StateVector[n]);
            }
            
            for (n=4;n<=7;n++) {
                StateVector[n] /= Math.sqrt(NormSquared[1]);
                //System.out.println("StateVector[" + n + "]" + StateVector[n]);
            }
            */
            
            
            for (n=0;n<=NormSquared.length-1;n++) {
                for (m=n*NumberOfRealValuesPerSite;m<=(n+1)*NumberOfRealValuesPerSite-1;m++) {
                    StateVector[m] /= Math.sqrt(NormSquared[n]);
                }
            }
            
            
            
            //System.out.println("oooooooooooooooooooooooooooooooooooo");
            
            return StateVector;
            
        }
        
        
        public double[] CalculateParticleDensities(double[] x) {
            /*
            int n;
            
            double[] rho   = new double[Nsites];
            
            double[] ReC0 = new double[Nsites];
            double[] ImC0 = new double[Nsites];
            double[] ReC1 = new double[Nsites];
            double[] ImC1 = new double[Nsites];
            
            double[] NormSquared = new double[Nsites];
            */
            
            /*
            for (n=0;n<=Nsites-1;n++) {
                ReC0[n] = x[4*n+0];
                ImC0[n] = x[4*n+1];
                
                ReC1[n] = x[4*n+2];
                ImC1[n] = x[4*n+3];
                
                NormSquared[n] =    ReC0[n]*ReC0[n] + ImC0[n]*ImC0[n] +
                                    ReC1[n]*ReC1[n] + ImC1[n]*ImC1[n]
                ;
                
                rho[n]=(ReC1[n]*ReC1[n]+ImC1[n]*ImC1[n]) /NormSquared[n];
            }
            */
            
            
            for (n=0;n<=Nsites-1;n++) {
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    ReC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m)  ];
                    ImC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m+1)];
                }
                
                NormSquared[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    NormSquared[n] +=   ReC[m][n]*ReC[m][n] + ImC[m][n]*ImC[m][n];
                }
                
                
                rho[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    rho[n] +=     (ReC[m][n]*ReC[m][n]+ImC[m][n]*ImC[m][n])
                                * (MinNumberOfParticles+m);
                }
                
                rho[n] /= NormSquared[n];
            }
            
            return rho.clone();
        }
        
        
        public double[] CalculateAbsoluteValueOfPsis(double[] x) {
            /*
            int n;
            
            double[] Repsi = new double[Nsites];
            double[] Impsi = new double[Nsites];
            
            double[] Abspsi = new double[Nsites];
            
            double[] ReC0 = new double[Nsites];
            double[] ImC0 = new double[Nsites];
            double[] ReC1 = new double[Nsites];
            double[] ImC1 = new double[Nsites];
            
            double[] NormSquared = new double[Nsites];
            */
            /*
            for (n=0;n<=Nsites-1;n++) {
                ReC0[n] = x[4*n+0];
                ImC0[n] = x[4*n+1];
                
                ReC1[n] = x[4*n+2];
                ImC1[n] = x[4*n+3];
                
                NormSquared[n] =    ReC0[n]*ReC0[n] + ImC0[n]*ImC0[n] +
                                    ReC1[n]*ReC1[n] + ImC1[n]*ImC1[n]
                ;
                
                Repsi[n] = (ReC0[n]*ReC1[n] + ImC0[n]*ImC1[n]) /NormSquared[n];
                Impsi[n] = (ReC0[n]*ImC1[n] - ImC0[n]*ReC1[n]) /NormSquared[n];
                
                Abspsi[n] = Repsi[n]*Repsi[n]+Impsi[n]*Impsi[n];
            }*/
            
            
            for (n=0;n<=Nsites-1;n++) {
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    ReC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m)  ];
                    ImC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m+1)];
                }
                
                NormSquared[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    NormSquared[n] +=   ReC[m][n]*ReC[m][n] + ImC[m][n]*ImC[m][n];
                }
                
                Repsi[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1-1;m++) {
                    Repsi[n] +=     (ReC[m][n]*ReC[m+1][n] + ImC[m][n]*ImC[m+1][n])
                                * Math.sqrt((MinNumberOfParticles+m)+1);
                }
                
                Repsi[n] /= NormSquared[n];
                
                
                Impsi[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1-1;m++) {
                    Impsi[n] +=     (ReC[m][n]*ImC[m+1][n] - ImC[m][n]*ReC[m+1][n])
                                * Math.sqrt((MinNumberOfParticles+m)+1);
                }
                
                Impsi[n] /= NormSquared[n];
                
                Abspsi[n] = Repsi[n]*Repsi[n]+Impsi[n]*Impsi[n];
                
            }
            
            
            return Abspsi;
        }
        
        
        public double[] NextDatapointForMinimization(double[] x) {
            
            double[] toreturn = new double[x.length];
            
            // coefficients for linear and quadratic term of quadratic approximation
            double a1,a2;
            
            double[] xMiddle;
            double[] xPlus;
            double[] xMinus;
            
            double xDifference12;
            double xDifference23;
            double xDifference31;
            
            double xDifferenceOfSquares12;
            double xDifferenceOfSquares23;
            double xDifferenceOfSquares31;
            
            double yMiddle;
            double yPlus;
            double yMinus;
            
            double Determinant;
            
            int n;
            
            //y = EnergyFunction.Value(x);
            
            for (n=0;n<=x.length-1;n++) {
            
                xMiddle = x.clone();
                yMiddle = energy(x);
                
                xPlus = x.clone();
                xPlus[n] += Math.pow(10.0, -2);//Math.pow(10.0, -10)
                yPlus = energy(xPlus);
                
                xMinus = x.clone();
                xMinus[n] -= Math.pow(10.0, -2);//Math.pow(10.0, -10)
                yMinus = energy(xMinus);
                
                xDifference12 = xMinus[n]  - xMiddle[n] ;
                xDifference23 = xMiddle[n] - xPlus[n]   ;
                xDifference31 = xPlus[n]   - xMinus[n]  ;
                
                xDifferenceOfSquares12 = xMinus[n]*xMinus[n]   - xMiddle[n]*xMiddle[n] ;
                xDifferenceOfSquares23 = xMiddle[n]*xMiddle[n] - xPlus[n]*xPlus[n]     ;
                xDifferenceOfSquares31 = xPlus[n]*xPlus[n]     - xMinus[n]*xMinus[n]   ;
                
                Determinant = (
                      xDifference12
                    * xDifference23
                    * xDifference31
                );
                
                if (
                    Determinant != 0.0
                ) {
                    a1 = (
                        (
                              yMinus  * xDifferenceOfSquares23
                            + yMiddle * xDifferenceOfSquares31
                            + yPlus   * xDifferenceOfSquares12
                        )
                        / Determinant
                    );
                    
                    a2 = (
                        (
                              yMinus  * (-1.0) * xDifference23
                            + yMiddle * (-1.0) * xDifference31
                            + yPlus   * (-1.0) * xDifference12
                        )
                        / Determinant
                    );
                    
                    /*
                    System.out.println( "-----------____________" );
                    
                    System.out.println( "x[" + n + "] = " + x[n] );
                    
                    System.out.println( "a1 n = " + n + ": " + a1 );
                    System.out.println( "a2 n = " + n + ": " + a2 );
                    
                    System.out.println( "yMinus = "  + yMinus  );
                    System.out.println( "yMiddle = " + yMiddle );
                    System.out.println( "yPlus = "   + yPlus );
                    
                    System.out.println( "yMiddle - yMinus = "  + (yMiddle - yMinus)  );
                    System.out.println( "yPlus - yMiddle = "   + (yPlus - yMiddle)   );
                    */
                    
                    if ( a2==0.0 ) { /* just linear approximation possible */
                        if (
                              ( ( xMiddle[n] - xMinus[n]  ) != 0.0 )
                            & ( ( xPlus[n]   - xMiddle[n] ) != 0.0 )
                        ) {
                            toreturn[n] = x[n] -
                                (
                                    (
                                        ( yMiddle - yMinus )/( xMiddle[n] - xMinus[n] )
                                    )
                                    +
                                    (
                                        ( yPlus   - yMiddle )/( xPlus[n]  - xMiddle[n] )
                                    )
                                    /2.0
                                )
                            ;
                        }
                        else {
                            toreturn[n] = x[n];
                        }
                    }
                    else if ( 0.0<=a2 ) {
                        if (
                              yMiddle < yMinus
                            & yMiddle < yPlus
                            & false
                        ) {
                            toreturn[n] = (-1.0)*a1/a2/2.0;
                            System.out.println( "\u001B[32m\u001B[1m-( a1/a2/2.0 ) = " + ( (-1.0)*a1/a2/2.0 ) + "\u001B[0m");
                        }
                        else {
                            if (
                                ( ( xMiddle[n] - xMinus[n]  ) != 0.0 )
                                & ( ( xPlus[n]   - xMiddle[n] ) != 0.0 )
                            ) {
                                toreturn[n] = x[n] -
                                    (
                                        (
                                            ( yMiddle - yMinus )/( xMiddle[n] - xMinus[n] )
                                        )
                                        +
                                        (
                                            ( yPlus   - yMiddle )/( xPlus[n]  - xMiddle[n] )
                                        )
                                        /2.0
                                    )
                                ;
                            }
                            else {
                                toreturn[n] = x[n];
                            }
                        }
                    }
                    else if ( a2<=0.0 ) {
                        if (
                              ( ( xMiddle[n] - xMinus[n]  ) != 0.0 )
                            & ( ( xPlus[n]   - xMiddle[n] ) != 0.0 )
                        ) {
                            toreturn[n] = x[n] -
                                (
                                    (
                                        ( yMiddle - yMinus )/( xMiddle[n] - xMinus[n] )
                                    )
                                    +
                                    (
                                        ( yPlus   - yMiddle )/( xPlus[n]  - xMiddle[n] )
                                    )
                                    /2.0
                                )
                            ;
                        }
                        else {
                            toreturn[n] = x[n];
                        }
                    }
                    
                }
                else {
                    toreturn[n] = x[n];
                }
                
                //System.out.println( "toreturn[" + n + "] = " + toreturn[n] );
                
            }
            
            return toreturn;
            
        }
        
        private double energy(double[] x) {
            
            double toreturn=0.0;
            
            /*
            int n;
            int m;
            
            double[] rho   = new double[Nsites];
            double[] Repsi = new double[Nsites];
            double[] Impsi = new double[Nsites];
            
            double[] ReC0 = new double[Nsites];
            double[] ImC0 = new double[Nsites];
            double[] ReC1 = new double[Nsites];
            double[] ImC1 = new double[Nsites];
            */
            
            double[][] InterSiteInteractionDistributionMatrix = 
            {
                { 0.0, 0.5, 0.5 },
                { 0.0, 0.0, 0.5 },
                { 0.0, 0.0, 0.0 }
            };
            /*
            double[] NormSquared = new double[Nsites];
            */
            /*
            for (n=0;n<=Nsites-1;n++) {
                ReC0[n] = x[4*n+0];
                ImC0[n] = x[4*n+1];
                
                ReC1[n] = x[4*n+2];
                ImC1[n] = x[4*n+3];
                
                NormSquared[n] =    ReC0[n]*ReC0[n] + ImC0[n]*ImC0[n] +
                                    ReC1[n]*ReC1[n] + ImC1[n]*ImC1[n]
                ;
                
                rho[n]=(ReC1[n]*ReC1[n]+ImC1[n]*ImC1[n]) /NormSquared[n];
                
                Repsi[n] = (ReC0[n]*ReC1[n] + ImC0[n]*ImC1[n]) /NormSquared[n];
                Impsi[n] = (ReC0[n]*ImC1[n] - ImC0[n]*ReC1[n]) /NormSquared[n];
                
            }
            */
            
            for (n=0;n<=Nsites-1;n++) {
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    ReC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m)  ];
                    ImC[m][n] = x[NumberOfRealValuesPerSite*n+(2*m+1)];
                }
                
                NormSquared[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    NormSquared[n] +=   ReC[m][n]*ReC[m][n] + ImC[m][n]*ImC[m][n];
                }
                
                
                rho[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    rho[n] +=     (ReC[m][n]*ReC[m][n]+ImC[m][n]*ImC[m][n])
                                * (MinNumberOfParticles+m);
                }
                
                rho[n] /= NormSquared[n];
                
                
                Repsi[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1-1;m++) {
                    Repsi[n] +=     (ReC[m][n]*ReC[m+1][n] + ImC[m][n]*ImC[m+1][n])
                                * Math.sqrt((MinNumberOfParticles+m)+1);
                }
                
                Repsi[n] /= NormSquared[n];
                
                
                Impsi[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1-1;m++) {
                    Impsi[n] +=     (ReC[m][n]*ImC[m+1][n] - ImC[m][n]*ReC[m+1][n])
                                * Math.sqrt((MinNumberOfParticles+m)+1);
                }
                
                Impsi[n] /= NormSquared[n];
                
                
                onsite[n]=0.0;
                
                for (m=0;m<=NumberOfCoordinatesPerSite-1;m++) {
                    onsite[n] +=     (ReC[m][n]*ReC[m][n]+ImC[m][n]*ImC[m][n])
                                * (MinNumberOfParticles+m) * ((MinNumberOfParticles+m)-1);
                }
                
                onsite[n] /= NormSquared[n];
                
                //System.out.println( "rho[" + n + "]" + rho[n] );
                
                
            }
            
            
            toreturn=0.0;
            
            for (n=0;n<=HoppingDistributionMatrixValue_RealPart.length-1;n++) {
            
                toreturn +=
                                (-1) * J
                            * HoppingDistributionMatrixValue_RealPart[n]
                            * 2.0
                            * (
                                      Repsi[HoppingDistributionMatrixIndex1[n]]*Repsi[HoppingDistributionMatrixIndex2[n]]
                                    + Impsi[HoppingDistributionMatrixIndex1[n]]*Impsi[HoppingDistributionMatrixIndex2[n]]
                                )
                ;
            
            }
            
            
            for (n=0;n<=HoppingDistributionMatrixValue_ImaginaryPart.length-1;n++) {
            
                toreturn +=
                                (-1) * J
                            * HoppingDistributionMatrixValue_ImaginaryPart[n]
                            * 2.0
                            * (
                                      Impsi[HoppingDistributionMatrixIndex1[n]]*Repsi[HoppingDistributionMatrixIndex2[n]]
                                    - Repsi[HoppingDistributionMatrixIndex1[n]]*Impsi[HoppingDistributionMatrixIndex2[n]]
                                )
                ;
            
            }
            
            
            for (n=0;n<=Nsites-1;n++) {
                for (m=0;m<=Nsites-1;m++) {
                    toreturn +=
                                  V/2.0
                                * InterSiteInteractionDistributionMatrix[n][m]
                                * 2.0
                                * (
                                        rho[n]*rho[m]
                                  );
                }
            }
            
            
            for (n=0;n<=Nsites-1;n++) {
                toreturn += 1.0/2.0 * onsite[n];
            }
            
            
            for (n=0;n<=Nsites-1;n++) {
                toreturn += (-1) * mu * rho[n];
            }
            
            
            toreturn /= (double)Nsites;
            
            return toreturn;
            
        }

    }

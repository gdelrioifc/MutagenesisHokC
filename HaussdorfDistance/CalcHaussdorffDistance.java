/*
Copyright (C) 2016 Victor Martinell, Gabriel Del Rio

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.algorithm.distance.DiscreteHausdorffDistance;

import java.util.Vector;
import java.util.Random;
import java.lang.Math;
import java.io.*;

public class CalcHaussdorffDistance 
{

	// Set initial temp
	private double temp;
	private double tempPercent = 0.15;
	private int searchStepsTotal;
	private int searchStepsCicle;

	// Cooling rate
	private double coolingRate;

	//factor de cambio para Phi y Psi
	private double initTheta;
	private double initTransX;
	private double initTransY;
	private double minTheta;
	private double minTransX;
	private double minTransY;
	private double cambioTheta;
	private double cambioTransX;
	private double cambioTransY;
	private double angleSteps = 0;

	private Geometry startGeometry;
	private Geometry targetGeometry;
	private String fileDir = "";

	private long initSeed;
	// 0 - RCC
	// 1 - RMSD
	private int energyType = 0;
	private boolean dualEnergy = false;
	private boolean rawAccept = false;

	String target_pdb = null;
	String current_pdb = null;
	//steps = -1/log((1/temp), 1+coolingRate)
	
	GeometryHandler gh;
	Random rdm;

	public CalcHaussdorffDistance(int searchStepsTotal, int searchStepsCicle, double anguloInicial, double anguloFinal, double transInicial, double transFinal, Geometry g1, Geometry g2, String fileDir, long initSeed){

		this.searchStepsTotal = searchStepsTotal;
		this.searchStepsCicle = searchStepsCicle;
		this.temp = 0;
		this.coolingRate = 0;

		this.initTheta = this.cambioTheta = anguloInicial;
		this.initTransX = this.cambioTransX = this.initTransY = this.cambioTransY = transInicial;
		this.minTheta = anguloFinal;
		this.minTransX = this.minTransY = transFinal;

		this.startGeometry = g1;
		this.targetGeometry = g2;
		this.fileDir = fileDir;
		this.initSeed = initSeed;
		rdm =  new Random(System.currentTimeMillis());
		gh = new GeometryHandler();
	}

	public void setEnergyType(int type){
		this.energyType = type;
	}

	public void setRawAccept(boolean rawAccept){
		this.rawAccept = rawAccept;
	}

	public double getTemp(){
		return this.temp;
	}

	public void setTemp(double temp){
		this.temp = temp;
	}

	public void setTempPercent(double tempPercent){
		this.tempPercent = tempPercent;
	}

	public void setCoolRate(double coolingRate){
		this.coolingRate = coolingRate;
	}

	//Calculate the acceptance probability
	public double acceptanceProbability(double energy, double newEnergy, double temperature) {
		if(rawAccept){
			return rawAcceptanceProbability(energy, newEnergy);
   	}
		// If the new solution is better, accept it
		if (newEnergy < energy) return 1.0;
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((energy - newEnergy) / temperature);
	}

	public double rawAcceptanceProbability(double energy, double newEnergy) {
		if (newEnergy < energy) return 1.0;
		return 0.0;
	}

	public double[] calcEnergy(Geometry g1, Geometry g2){
		double energy[] = {0.0,0.0};
		energy[0] = 100* DiscreteHausdorffDistance.distance(g1,g2);
		return energy;
	}

	public Geometry alterConformation(Geometry inGeo){
		Random rdm =  new Random(System.currentTimeMillis());
		double dTheta = cambioTheta * 2*(rdm.nextDouble() - 0.5);
		double dX = cambioTransX * 2*(rdm.nextDouble() - 0.5);
		double dY = cambioTransY * 2*(rdm.nextDouble() - 0.5);
		Geometry alterGeo;
		alterGeo = gh.rotate(inGeo, dTheta);
		alterGeo = gh.translate(alterGeo, dX, dY);
		return alterGeo;
	}

	//
	//Random de 0 error
	//al dividir salen 0 donde no deberían
	//todo mayor que 0!!!!!

	public double angleCooling(double time){
		return initTheta*Math.exp(-time/(-1.0/Math.log(minTheta/initTheta)));
	}

	public double transCooling(double time){
		return initTransX*Math.exp(-time/(-1.0/Math.log(minTransX/initTransX)));
	}

	public double stdCooling(double temp){
		return temp * (1-coolingRate);
	}

	public double linearCooling(double temp){
		return temp - coolingRate;
	}

	public double slowCooling(double temp){
		return temp/Math.log(temp);
	}

	//verbos 0 - none
	//			 1 - just solutions
	//			 2 - all
	//			 3 - temp calc
	
	//public double calcInitialTemp(int verbos){
	//	int count = 0;
	//	int searchStepsCicle = 100;
	//	double sum = 0;
	//	double currentEnergy[]={0.0,0.0};
	//	while(count < searchStepsCicle){
	//		Geometry  = alterConformationAll(alterStruc);

	//		sum += currentEnergy[energyType];
	//		if(verbos > 2){
	//			System.out.println("energy " + currentEnergy[energyType] + " sum " + sum);
	//		}
	//		count++;
	//	}
	//	double prom = sum/count;
	//	if(verbos > 1){
	//		System.out.println(" promE " + prom);
	//	}
	//	double temp = prom;//2*max - prom;
	//	this.temp = temp * tempPercent;
	//	return temp;
	//}

	public void calcCoolRateStd(int verbos, int searchStepsTotal, int searchStepsCicle){
		double steps = searchStepsTotal/(double)searchStepsCicle;
		this.coolingRate = 1-(Math.pow((1/this.temp),(1.0/steps)));
	}
	
	public String initialize(int verbos){
		// Load Protein Sequence
		String outString = null;
		
		// Declaration of variables to store energy values
		// TEMP
	
		//if(temp == 0){
		//	calcInitialTemp(verbos, targetGeometry);
		//}
		if(coolingRate == 0){
			calcCoolRateStd(verbos, searchStepsTotal, searchStepsCicle);
		}

		if (verbos > 0){
			if (verbos > 1){
				System.out.print("\n");
			}
			return outString;
		}
		return null;
	}

	public double run(int verbos){
		//Mesure time
		long elapsedTime = 0;
		if (verbos > 0){
			elapsedTime = System.nanoTime();
		}
		
		Geometry struc_fit = startGeometry;
		Geometry struc_model = null;
		double currentEnergy[] = {0.0,0.0};
		double neighbourEnergy[] = {0.0,0.0};
		double distance_ini[] = {0.0,0.0};
		double best;
		
		// Initialize intial solution
		distance_ini = calcEnergy(startGeometry,targetGeometry);
		
		if (verbos > 0){
			System.out.println("Initial solution distance: " + distance_ini[0] + " " + distance_ini[1]);
		}
		
		// Set as current best
		best = distance_ini[energyType];
		    
		// Create new neighbour 3d model
		currentEnergy = distance_ini;

		// Loop until system has cooled
		for (;temp > 1;temp = stdCooling(temp)){
			// Search steps
			cambioTheta = angleCooling(angleSteps);
			cambioTransX = transCooling(angleSteps);
			cambioTransY = cambioTransX;
			/// REVISARAR
			angleSteps += (double)searchStepsCicle/searchStepsTotal;
			for(int step = 0; step < searchStepsCicle; step++){
				// Get a random conformation for this new neighbor
				struc_model = alterConformation(struc_fit);
				// Get energy of solution
				neighbourEnergy = calcEnergy(struc_model,targetGeometry);
			
				// Decide if we should accept the neighbour
				if (acceptanceProbability(currentEnergy[energyType], neighbourEnergy[energyType], temp) > rdm.nextDouble()) {
					struc_fit = struc_model;
					currentEnergy = neighbourEnergy;
					if (verbos > 0){
						if (verbos > 1){
							long time = System.nanoTime() - elapsedTime;
							System.out.println(temp + "> " + currentEnergy[energyType] + " time: " + time/3600000000000.0);
						}else{
							System.out.println(temp + "> " + currentEnergy[energyType]);
						}
					}
				}
			
				// Keep track of the best solution found
				if (neighbourEnergy[energyType] < best){
					best = neighbourEnergy[energyType];
					if(best == 0){
						temp = 0;
					}
				}
				if (verbos > 1){
					//if(energyType == 0 || dualEnergy){
					//	for(int i : modelRCC){System.out.print(i + " ");}
					//	System.out.print("\n");
					//}
					if(dualEnergy){
						System.out.println(neighbourEnergy[0] +"\t\t"+ neighbourEnergy[1]);
					}else{
						System.out.println(neighbourEnergy[energyType] +"\t\t"+ cambioTheta);
					}
				}
			}
		}
		
		if (verbos > 0){
			System.out.println("Final solution distance: " + best);//.getDistance());
			elapsedTime = System.nanoTime() - elapsedTime;
			System.out.println("Total execution time: " + elapsedTime/3600000000000.0);
		}
		return best;
  }

	public static void main(String[] args){
		if(args.length<2 || args.length>3){
     System.err.println("Gnomovision version 69, Copyright (C) 2016 Victor Martinell and Gabriel Del Rio");
     System.err.println("Gnomovision comes with ABSOLUTELY NO WARRANTY.");
     System.err.println("Please read the conditions at http://www.gnu.org/licenses/gpl-2.0.html.\n\n");

			System.err.println("Usage:\njava CalcHaussdorffDistance <sphCoorRef> <sphCoor> <outputDirectory>");
 			System.err.println("<sphCoorRef> and <sphCoor> files generated with CalcSphCoordPDB.");
			System.err.println("<outputDirectory> directory name where to save the output; if no <outputDirectory> is provided, out/ will be used.");
			System.err.println("Requires external package for compiling and executing; please download it from: http://www.vividsolutions.com/jts/JTSHome.html");
			System.err.println("The program will compare the two sphere coordinates provided in <sphCoorRef> and <sphCoor> to compute the minimum");
			System.err.println("Haussdorff distance, by using an simulated annealing algorithm.");
		}else{
		int searchStepsTotal = 300000;
		int searchStepsCicle = 20;

		double anguloInicial = 0.0005;
		double anguloFinal = 0.00001;
		double transInicial = 0.00005;
		double transFinal = 0.000001;

		String fileDir = "out/";
		long initSeed = (long)(1000*Math.random());
		if(args.length>2){
			fileDir = args[2];
		}
		GeometryHandler f = new GeometryHandler();
		Geometry g1 = f.toGeometry(f.readCoordinate(args[0]));
		Geometry g2 = f.toGeometry(f.readCoordinate(args[1]));
		CalcHaussdorffDistance simA = new CalcHaussdorffDistance(searchStepsTotal,searchStepsCicle,anguloInicial,anguloFinal,transInicial,transFinal,g1,g2,fileDir,initSeed);
		////SimulatedAnnealing3DProtFromRCC simA = new SimulatedAnnealing3DProtFromRCC(args[0], args[1]);
		//simA.setRawAccept(true);
		simA.setTemp(2.5);
		simA.initialize(2);
		simA.run(2);
		}
	}
}

/*
Copyright (C) 2016 Gabriel Del Rio

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


import java.io.*;
import java.util.*;

/**
 *  Description of the Class
 *
 *@author     gdelrio
 *@created    July 25, 2003
 */
public class doGraphFromPDB {
	static LinkedList atoms = new LinkedList();
	static Hashtable atoms_aa_num = new Hashtable();
	static Hashtable atoms_ht = new Hashtable();
	static Hashtable aa_num_avgs = new Hashtable();
	static Vector p_charged_aa = new Vector();
	static Vector n_charged_aa = new Vector();
	static Vector hydrophobs_aa = new Vector();


	/**
	 *  Description of the Method
	 */
	public static void initializeAA() {
		p_charged_aa.add("ASP");
		p_charged_aa.add("GLU");
		p_charged_aa.add("ASN");
		p_charged_aa.add("GLN");
		p_charged_aa.add("TYR");
		n_charged_aa.add("ARG");
		n_charged_aa.add("LYS");
		n_charged_aa.add("ASN");
		n_charged_aa.add("GLN");
		n_charged_aa.add("PHE");
		n_charged_aa.add("HIS");
		hydrophobs_aa.add("ALA");
		hydrophobs_aa.add("LEU");
		hydrophobs_aa.add("ILE");
		hydrophobs_aa.add("VAL");
		hydrophobs_aa.add("TRP");
		hydrophobs_aa.add("SER");
		hydrophobs_aa.add("THR");
		hydrophobs_aa.add("CYS");
		hydrophobs_aa.add("MET");
	}


	/**
	 *  Gets the molecule attribute of the doGraphFromPDB class
	 *
	 *@param  file                       Description of the Parameter
	 *@param  chain                      Description of the Parameter
	 *@exception  IOException            Description of the Exception
	 *@exception  FileNotFoundException  Description of the Exception
	 */
	public static void getMolecule(String file, String chain) throws IOException, FileNotFoundException {
		String line = "";
		String resName = "";
		String chainID = "";
		String resSeq = "";
		String serial = "";
		String xS = "";
		String yS = "";
		String zS = "";
		String aa = "";
		String aa_old = "";
		double avgX = 0.0;
		double avgY = 0.0;
		double avgZ = 0.0;
		BufferedReader infile = new BufferedReader(new FileReader(file));
		Double[] coords = null;
		Double[] avg_coords = null;
		int i = 0;

		while ((line = infile.readLine()) != null) {
			if (line.startsWith("ATOM")) {
				chainID = line.substring(21, 22);
				if (chainID.equalsIgnoreCase(chain) || chain.equalsIgnoreCase("none") || chain.equalsIgnoreCase("all")) {
					resName = line.substring(17, 20);
					resSeq = line.substring(22, 26);
					serial = line.substring(6, 11);
					xS = line.substring(30, 38);
					yS = line.substring(38, 46);
					zS = line.substring(46, 54);

					coords = new Double[3];
					coords[0] = new Double(Double.parseDouble(xS));
					coords[1] = new Double(Double.parseDouble(yS));
					coords[2] = new Double(Double.parseDouble(zS));

					if (chain.equalsIgnoreCase("all")) {
						aa = resName.trim() + resSeq.trim() + "_" + chainID.toUpperCase();
					} else {
						aa = resName.trim() + resSeq.trim();
					}

					if (!atoms.contains(serial)) {
						atoms.add(serial);
					}

					if (!atoms_ht.containsKey(serial)) {
						atoms_ht.put(serial, coords);
					}

					if (!atoms_aa_num.contains(serial)) {
						atoms_aa_num.put(serial, aa);
					}

					if (aa_num_avgs.containsKey(aa)) {
						avgX = avgX + Double.parseDouble(xS);
						avgY = avgY + Double.parseDouble(yS);
						avgZ = avgZ + Double.parseDouble(zS);
						i++;
						aa_old = aa;
					} else {
						if (!aa_num_avgs.isEmpty()) {
							avg_coords = new Double[3];
							avg_coords[0] = new Double(avgX / i);
							avg_coords[1] = new Double(avgY / i);
							avg_coords[2] = new Double(avgZ / i);

							aa_num_avgs.remove(aa_old);
							aa_num_avgs.put(aa_old, avg_coords);

						}
						avgX = Double.parseDouble(xS);
						avgY = Double.parseDouble(yS);
						avgZ = Double.parseDouble(zS);
						aa_num_avgs.put(aa, new Double[3]);
						i = 1;
					}
				}
			}
		}

		avg_coords = new Double[3];
		avg_coords[0] = new Double(avgX / i);
		avg_coords[1] = new Double(avgY / i);
		avg_coords[2] = new Double(avgZ / i);

		aa_num_avgs.remove(aa_old);
		aa_num_avgs.put(aa_old, avg_coords);

		infile.close();

	}
	// end of getMolecule method


	/**
	 *  Description of the Method
	 *
	 *@param  aa1    Description of the Parameter
	 *@param  aa2    Description of the Parameter
	 *@param  pairs  Description of the Parameter
	 *@return        Description of the Return Value
	 */
	public static boolean complement(String aa1, String aa2, String pairs) {
		boolean result = false;

		if (pairs.equalsIgnoreCase("all_complementary")) {
			if ((hydrophobs_aa.contains(aa1) && hydrophobs_aa.contains(aa2)) || (p_charged_aa.contains(aa1) && n_charged_aa.contains(aa2)) || (p_charged_aa.contains(aa2) && n_charged_aa.contains(aa1))) {
				result = true;
			}
		} else if (pairs.equalsIgnoreCase("h_complementary")) {

			if (hydrophobs_aa.contains(aa1) && hydrophobs_aa.contains(aa2)) {
				result = true;
			}
		} else if (pairs.equalsIgnoreCase("ch_complementary")) {
			if ((p_charged_aa.contains(aa1) && n_charged_aa.contains(aa2)) || (p_charged_aa.contains(aa2) && n_charged_aa.contains(aa1))) {
				result = true;
			}
		} else {
			return false;
		}
		return result;
	}
	// end of complement method


	/**
	 *  Description of the Method
	 *
	 *@param  min                        Description of the Parameter
	 *@param  max                        Description of the Parameter
	 *@param  distance_criterium         Description of the Parameter
	 *@return                            Description of the Return Value
	 *@exception  IOException            Description of the Exception
	 *@exception  FileNotFoundException  Description of the Exception
	 */
	public static String doDistance(double min, double max, String distance_criterium) throws IOException, FileNotFoundException {
		String num1 = "";
		String num2 = "";
		String token = "";
		String root_file_name = "distance_file_";
		String file_name = "";
		String aa1 = "";
		String aa2 = "";
		String aa1_old = "";
		String aa2_old = "";
		int i = 0;
		int j = 0;
		boolean file_test = false;
		double x1 = 0.0;
		double y1 = 0.0;
		double z1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;
		double z2 = 0.0;
		double xx = 0.0;
		double yy = 0.0;
		double zz = 0.0;
		double distance = 0.0;
		Double[] coords = null;
		Double[] avg_coords = null;
		File file = null;
		LinkedList neighbors = new LinkedList();

		while (file_test == false) {
			file_name = root_file_name + String.valueOf(i);
			file = new File(file_name);

			if (file.exists()) {
				i++;
			}
			if (!file.exists()) {
				file_test = true;
			}
		}

		BufferedWriter outfile = new BufferedWriter(new FileWriter(file_name));

		if (distance_criterium.equalsIgnoreCase("once")) {
			ListIterator l1 = atoms.listIterator(0);
			while (l1.hasNext()) {
				num1 = (String) l1.next();
				aa1 = (String) atoms_aa_num.get(num1);
				if (aa1_old.equals("")) {
					aa1_old = aa1;
				}

				coords = (Double[]) atoms_ht.get(num1);
				x1 = coords[0].doubleValue();
				y1 = coords[1].doubleValue();
				z1 = coords[2].doubleValue();

				if (!aa1.equals(aa1_old)) {
					if (neighbors.size() > 0) {
						for (j = 0; j < neighbors.size(); j++) {
							outfile.write(aa1_old + "\t" + (String) neighbors.get(j) + "\n");
						}
					}
					neighbors.clear();
					aa1_old = aa1;
				}
				ListIterator l2 = atoms.listIterator(l1.nextIndex());
				while (l2.hasNext()) {
					num2 = (String) l2.next();
					aa2 = (String) atoms_aa_num.get(num2);

					coords = (Double[]) atoms_ht.get(num2);
					x2 = coords[0].doubleValue();
					y2 = coords[1].doubleValue();
					z2 = coords[2].doubleValue();

					xx = Math.pow((x1 - x2), 2.0);
					yy = Math.pow((y1 - y2), 2.0);
					zz = Math.pow((z1 - z2), 2.0);
					distance = Math.sqrt(xx + yy + zz);
					if (!aa1.equalsIgnoreCase(aa2)) {
						if (distance >= min && distance <= max) {
							if (!neighbors.contains(aa2)) {
								neighbors.add(aa2);
							}
						}
					}
				}
			}
			for (j = 0; j < neighbors.size(); j++) {
				outfile.write(aa1 + "\t" + (String) neighbors.get(j) + "\n");
			}
		}

		if (distance_criterium.equalsIgnoreCase("average")) {
			for (Enumeration e1 = aa_num_avgs.keys(); e1.hasMoreElements(); ) {
				aa1 = (String) e1.nextElement();
				avg_coords = (Double[]) aa_num_avgs.get(aa1);
				x1 = avg_coords[0].doubleValue();
				y1 = avg_coords[1].doubleValue();
				z1 = avg_coords[2].doubleValue();

				for (Enumeration e2 = aa_num_avgs.keys(); e2.hasMoreElements(); ) {
					aa2 = (String) e2.nextElement();
					if (!aa1.equalsIgnoreCase(aa2)) {
						avg_coords = (Double[]) aa_num_avgs.get(aa2);
						x2 = avg_coords[0].doubleValue();
						y2 = avg_coords[1].doubleValue();
						z2 = avg_coords[2].doubleValue();

						xx = Math.pow((x1 - x2), 2.0);
						yy = Math.pow((y1 - y2), 2.0);
						zz = Math.pow((z1 - z2), 2.0);
						distance = Math.sqrt(xx + yy + zz);

						if (distance >= min && distance <= max) {
							outfile.write(aa1 + "\t" + aa2 + "\n");
						}
					}
				}
			}
		}

		outfile.flush();
		outfile.close();
		return file_name;
	}
	// end of doDistance method


	/**
	 *  Description of the Method
	 *
	 *@param  file                       Description of the Parameter
	 *@param  pairs                      Description of the Parameter
	 *@return                            Description of the Return Value
	 *@exception  IOException            Description of the Exception
	 *@exception  FileNotFoundException  Description of the Exception
	 */
	public static boolean printResults(String file, String pairs, String outname) throws IOException, FileNotFoundException {
		boolean result = false;
		BufferedReader infile = new BufferedReader(new FileReader(file));
		BufferedWriter outfile = new BufferedWriter(new FileWriter(outname));
		File to_be_deleted = new File(file);
		String line = "";
		String res1 = "";
		String res2 = "";
		String aa1 = "";
		String aa2 = "";
		StringTokenizer st = null;

		while ((line = infile.readLine()) != null) {
			st = new StringTokenizer(line, "\t");
			if (st.countTokens() == 2) {
				res1 = (String) st.nextElement();
				res2 = (String) st.nextElement();

				if (pairs.equalsIgnoreCase("all")) {
					outfile.write(res1 + "\t" + res2 + "\n");
				} else {
					aa1 = res1.substring(0, 3);
					aa2 = res2.substring(0, 3);

					if (complement(aa1, aa2, pairs)) {
						outfile.write(line+"\n");
					}
				}
			}
			while (st.hasMoreElements()) {
				st.nextElement();
			}
		}
		infile.close();
		outfile.flush();
		outfile.close();

		if (to_be_deleted.delete()) result = true;

		return result;
	}
	// end of printResults method


	/**
	 *  The main program for the doGraphFromPDB class
	 *
	 *@param  args                       The command line arguments
	 *@exception  IOException            Description of the Exception
	 *@exception  FileNotFoundException  Description of the Exception
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException {
		if (args.length == 7) {
			String pdb = args[0];
			String chain = args[1];
			String pairs = args[4];
			String distance_criterium = args[5];
			String outfile = args[6];
			String num1 = "";
			String aa1 = "";
			String num2 = "";
			String aa2 = "";
			String file_name = "";
			double min = Double.parseDouble(args[2]);
			double max = Double.parseDouble(args[3]);

			System.err.print("Initizalizing variable...");
			initializeAA();
			System.err.println("Done.");
			System.err.print("Loading " + pdb + " ...");
			getMolecule(pdb, chain);
			System.err.println("Done");
			System.err.print("Calculating interatomic distances...");
			file_name = doDistance(min, max, distance_criterium);
			System.err.println("Done");
			System.err.print("Writing and formating results...");
			while (printResults(file_name, pairs, outfile) == false) {
				;
			}
			System.err.println("Done.");
		} else {
			System.err.println("Usage:\njava doGraphFromPDB <pdb> <chain> <min> <max> <pairs> <distance_criterium> <output_filename>");
			System.err.println("<pdb> pdb file name");
			System.err.println("<chain> chain name in <pdb> to be used to trace the graph. If no chain is defined");
			System.err.println("in the <pdb>, you must make <chain>=\"none\". If you want to consider all, <chain>=\"all\".");
			System.err.println("<min> float specifying minimum distance to find neighbors");
			System.err.println("<max> float specifying maximum distance to find neighbors");
			System.err.println("<distance_criterium> is either \"average\" or \"once\".");
			System.err.println("<pairs> is either \"all\", \"all_complementary\", \"h_complementary\" or \"ch_complementary\".");
			System.err.println("The program will trace a graph of amino acid residues neighbors");
			System.err.println("from <pdb> using <min> and <max> as a cutoff criteria. <pairs> specifies which residues");
			System.err.println("within <min> and <max> will be paired:");
			System.err.println("\"all\": includes all residues.");
			System.err.println("\"all_complementary\": includes complementary charged (positive with negative) and hydrophobic residues");
			System.err.println("\"h_complementary\": includes complementary hydrophobic residues");
			System.err.println("\"ch_complementary\": includes complementary charged residues");
			System.err.println("<distance_criterium> determines whether two residues are neighbors. If \"once\",");
			System.err.println("two residues are neighbors if at least one atom in each residue is within <min><max> range.");
			System.err.println("If \"average\", two residues are neighbors if the average distance between the two");
			System.err.println("is within the <min>max> range.");
			System.err.println("<output_filename> the name of a file where to write the resulting graph.");
			System.err.println("The output will be a two column tab-delimited.");
			System.err.println("An additional output will be sent to the error standard output reporting the estimated distances.");
			System.err.println("");
		}
	}
}


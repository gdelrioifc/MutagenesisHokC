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

public class calcSizeOrderGraph
{
 static LinkedList vertices = new LinkedList();
 static LinkedList edges = new LinkedList();

 public static void main(String[] args) throws Exception
  {
   if(args.length==1)
    {
     String line="", v1="", v2="", file=args[0];
     StringTokenizer st = null;
     LinkedList edge = null;
     BufferedReader infile = new BufferedReader(new FileReader(file));

     while((line=infile.readLine())!=null)
      {
       st = new StringTokenizer(line, "\t");
       if(st.countTokens()>=2)
	{
         v1 = (String)st.nextElement();
         v2 = (String)st.nextElement();

	 if(!vertices.contains(v1)) vertices.add(v1);
	 if(!vertices.contains(v2)) vertices.add(v2);

	 edge = new LinkedList();
	 edge.add(v1); edge.add(v2);
	 Collections.sort(edge);
	 if(!edges.contains(edge)) edges.add(edge);
	}
       while(st.hasMoreElements()) st.nextElement();
      }
     infile.close();

     System.out.println("Size= "+String.valueOf(vertices.size())+"\nOrder= "+String.valueOf(edges.size()));
    }
   else
    {
     System.err.println("Gnomovision version 69, Copyright (C) 2016 Gabriel Del Rio");
     System.err.println("Gnomovision comes with ABSOLUTELY NO WARRANTY."); 
     System.err.println("Please read the conditions at http://www.gnu.org/licenses/gpl-2.0.html.\n\n");

     System.err.println("Usage:\njava calcSizeOrderGraph <file>");
     System.err.println("<file> file with two columns tab-delimited, the graph.");
     System.err.println("The program will count the number of unique vertices and edges.");
     System.err.println("The output will be sent to the standard output.");
     System.err.println("");
    }
  }
}

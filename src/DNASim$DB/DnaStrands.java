package DNASim$DB;

import java.util.Arrays;

import TaiGameCore.GameDataBase;

public class DnaStrands extends GameDataBase{
	public DnaStrands(String hash) {
		super(hash);
	}
	public static void main(String[] args){
		//Writes a DNA Strand.
		DnaStrands made = new DnaStrands("");
		made.strands = new DnaStrand[6];
		made.domains = new DnaDomain[14];
		made.domainComplements = new int[14];
		int i = 0; 
		for(String domain : "CatcgcacctaagcatactC CctagcacatatagccgcaC CccgaaacttgaatttctcC CttactcgcgtaaactaacC CcttttcaaatgtactacaC CtcaagatactccattaagC CttcagcaatcaacagcaaC".split(" ")){
			domain = domain.replace("\\s","");
			DnaDomain newDomain = new DnaDomain("");
			DnaDomain newDomainC = new DnaDomain("");
			
			int numChars = domain.length();
			int charInd = 0;
			newDomain.data = new int[numChars];
			newDomainC.data = new int[numChars];
			for(char q : domain.toCharArray()){
				newDomain.data[charInd] = DnaDomain.readNA(q, false);
				newDomainC.data[charInd] = DnaDomain.readNA(q, true);
				charInd++;
			}
			made.domainComplements[i*2] = i*2+1;
			made.domainComplements[i*2+1] = i*2;
			made.domains[i*2] = newDomain;
			made.domains[i*2+1] = newDomainC;
			i++;
		}
		final int COMP = 1;
		DnaStrand strand;
		//Just go ahead and make strands too.
		DnaStrand OP = new DnaStrand("");strand = OP;
		strand.domains = new int[]{1*2,2*2};
		strand.name = "OP";
		makeCopies(strand,1);
		
		DnaStrand SUBST = new DnaStrand("");strand = SUBST;
		strand.domains = new int[]{2*2+COMP,3*2+COMP,4*2+COMP,5*2+COMP};
		strand.name = "SUBST";
		makeCopies(strand,1);
		
		DnaStrand SIDE = new DnaStrand("");strand = SIDE;
		strand.domains = new int[]{6*2,3*2,4*2};
		strand.name = "SIDE";
		makeCopies(strand,1);
		
		DnaStrand FUEL = new DnaStrand("");strand = FUEL;
		strand.domains = new int[]{2*2,3*2,4*2};
		strand.name = "FUEL";
		makeCopies(strand,1);
		
		DnaStrand CAT = new DnaStrand("");strand = CAT;
		strand.domains = new int[]{4*2,5*2};
		strand.name = "CAT";
		makeCopies(strand,1);
		
		made.strands = new DnaStrand[]{
				OP,SUBST,SIDE,FUEL,CAT
		};
		System.out.print(made.hashToString());
	}
	private static void makeCopies(DnaStrand stra, int i) {
		stra.positions = new float[stra.numDimensions * i * stra.domains.length];
		stra.velocities = new float[stra.numDimensions * i * stra.domains.length];
		stra.bindersInd = new int[stra.numDimensions * i * stra.domains.length];
		stra.bindersStrand = new int[stra.numDimensions * i * stra.domains.length];
		Arrays.fill(stra.bindersInd,-1);
		Arrays.fill(stra.bindersStrand,-1);
		for(int y = 0; y < stra.positions.length; y+=stra.numDimensions*stra.domains.length){
			for(int k = 0; k < 3; k++){
				float got = (float)Math.random()-.5f;
				for(int q = 0; q < stra.domains.length; q++){
					stra.positions[y+k+q*stra.numDimensions] = (float) (got + Math.random()*.05f);
				}
			}
		}
		stra.numCopies = i;
	}
	public static DnaStrands makeFile() {
		return new DnaStrands("A3O1GJqJqLAJIJAJIJAJC3E3Ev0PAJAJAJAJAJAJAJAJAJAJAJAJCJAJAJAJC12vODK1EDOvurCLuBOFsNOrGF0vsD6BCvsJIPCBO5EJq1GBqNIJEDwJAturCLEtErCL4DIJALAJGrwNCxCN07sP6tC1K363uNArCvCxK1EB4LKtCrwvu7w1AB05AJsPu52xsxwt6DK3sr6x0xwt27sx6Nwv47Cv434HA76tMHwDwLwx6J6xCxwP6xCxw56xCxwL6xCxw76xuP4H07Ox616xCxw36xCx23wx0vOPO10xEv4J4xGFqD034r2xC14FwLwF6763Mrsv6D27sv07sD4r2xM14F63MrsxOLKH472Nwt6vqvMx2F036HwDOvO5Krq3M1KN4v6BGNwPG7GNwHs7GJwDMHsDOtG72FsvKtMr2NsvMP6x4BsN23Ers3OHAr2vI7utwvMB6PIH6LI12HGNKFwxIvMNwPwDOxIrsD45CNs561sF4B2D0P6tItsPA7KHqF61KrwPwDux0JI5w50N2v0Nsv2x0r21qN6xMtEt4FK1wFw3q525wx6ru7E3uF6rGHGDw5CDOLMHE1MvK3K14xKJ6vu5EDwPGDwvCPKvuJ07OxwvO7KF6NqDqHO7K7w3K1MDAD07CD0t6HMHCHqv6JqvOPqDExwDE341EN256JsN21ANGJqDIrGxqBuH2FsF0FKHqtwHEv4twxONuF2xGJKr6HuF2N4H6P2B0H6tM14H2Nux6P0DKtK1wruFCHC7I5qN67ANsv2rI1wHCN6vKH21GPO7w72D6PuF4xGtIt4xqtCP4vO36tsxIL6x0tONuF4L25O7KNwxst6JAvKF2x4JCHwxKv27GNIPMHs3ux6D4t4LCxADArOLGr4H2v0t0BEPw3uJKvsPADuxAx63EJGv2JuF21CBEHMD0705At6rsN23A5K7ELKx6PwrAJOJK76v2rMxMNG163Ir6xOJKx67MN2xKBqv27qN6xONqHONODu5O1uF4BuHKPGN6tuBK5APqF6tO5uFKrIDK5w7OL4vsFGxqxq7q5M1M7Otq36JuBO76tMx4DKJ4vsvGxAr2x2BCD47ENIxGLC7OvGP6BA76BqxsPIBsxsxM7ELwxIFwxIv4rqvKxALE56DAB0tq1IBwHCDO7s1sv2tK123K1AvMLMF21wPAJw1K5AxEBw5sL2L670NqtOxutwx01A36t6HO3IFwFuPMtGP4v6BE347CJq7GJA54vK1AFuJ4x4P670NAtGxqF4vC50HwHw54DIrwrIvML636tCD0HKNA34DM1MxMtADAHuF2x4LurCJG3G1Cv4Ds7A7wt07O7MxEvKB4xANs54BIDE5s3wPE3wDIv41AtC3wJMrIDCt6BwrEFOH652765uB4vu3KJ0xILqrGDw7OJCxGLODO7IJON0xu5APCHw5I56Jwx2J0N6vK3O7KB47KLAt4DqDE7sP21A50J0tEr47uNA3ItOrOxA3Mt6HqDuPwD0J6rArsHM7q5OB4NMLADuNKvAF6HOvAF6Pu5INwLEJw56r4v6FuHEF07IJG76PstE5urqNEF21OxwJGr0H65276rq1EvKvOHs1I7MtEDA145KBwx0J6xqJI16vCLOLMJwBu3EDwJArC7GNAPG361qt63KLu1KDGJCBKL0v0xKNE7072JGtwBuPs3u52DCP6NCv2vOP2PqF2D47MLu101Krw1ErOH65276rqJ6vIx2tErOxMP4tEDMD4tI1CDKr4D0H6BuNwDw7uFw5utI1K56FG1MLGNMxq3MtKtABMxEBw5O7Atwv2LEtMLKP2tI1OHsDIBw32LC7wF2LC76J03O10t2BO5sPOrC14vOJK5wxCDu7sDut4vO3IN6rqF0N6F2B0HGFODGvqL4PqJsvM1EJ2vAvux2FqxKDAF4x4HsvI7OLG5KB6r0H616JKvAL21CHE1ONqru1qtOr4NuFEFCxsD0DKFCPC1EDONKP6tO1GF01Irw7OJ4tO3G1G3A1I565GxA7EL6vO36tGPCNOxCvqN472LGLI74NI76tMJ0FqvEx0FE5ItqPM34Lw7EP274JqtstANG7qH01E5C1sDqHGt2N6Fut4twHOBIDOxqPKHIL6PGru3GFwBwBALOxutwx01A7676HCNqrG303srGF2vsxErCx052DMrGDGvKPArCx0L4D2vsJsPuFIrGPK7wPq56DqFwJKDONKJ0xqLu16H2vs7u5EDuFED0HO3qv6xAvCvCrsD4Hsx2PGLqr6x6Lw3KHqDsv2NC3w1uvAB4F4rAF2vKJCvsLqNuDuvqLuN05GDuLuD6t6vCJst01sDO7OrGxCvsB2D4NCt4HwPO3K7O7q5012HwrEF2H65276PGL2vwHsDAr2PAN632HO3Ir63OrGx0rE50NMLuP0DAFK1ELOBEDO7O3qvAPIxIr0HC32BM3wJqDOHuJ4tG5OBKH0t6NutEPE5MN2tEBKH6NwNqruP47q1wL6J6JsNux2F6vKDI12x67uNqDIFqNA76HABO1qruv0BMJwDKNuFst2L6x23wJKvuL4tG5MvMDENI54tM5G3I561ODOJuFqrADu70PMDKrOxstwx01AL676HMHsvu30Bs5K5wxwFG72NOLE54NOt45w367ArKHED0BIP01q5w5w34vO10tKvw1A323C7uBuvsLKv6H2vu16HwJMFsDqJMxEHw76HM16DKJMtGHOtAtuxwNqLOHGt4Bs1w7u5EDGLAD6tMvKJI1K1Exsx6H6D45EJqvq1u1MtOvIBMFI1CxKxq7EPqN2LIr61GLqFErw5q3AxsF2LIxuP2LIxMr4r630PCtA7ONItONOr63C5OtCJMFuNAN6v0D2LGFKD234303OxAFq5G1OxwNIBExKLM5I5s5CL2L25G50BO5w3O3K3OLOvwvuJM16rCH2D6L6NsF67sLOxwJODI1G3K1GPEHKJMPED0xOP0DIBqDKJwHErKBI5A7ILENEN4PsLq5AvGx6BA34LKPG5wF47uHwxGFG1uD0FIr0v2HCrCvItGH0525qJOvsBsF2FKFGFC54JI52PG50LMDKJKvuxqFGxIL27sDO725w7MrIHqLAvs7016xMt47IvCrIHqLE5IDE5OF0t4twBIxKxsHODKDEJ0F05GJM3OtMPEF4tuBKPuPArCxC3OPuLED6xsDExuBqx2Fu5KB2tE7Mvu541qJ2NA50JwFEx47C5wvuL0xIL0L2tI3K1ArOLI5wN6NOtGFMPOJuNM745EPutOL4t6vKLqP0Lq7IDE7AvC5ONsLuDAB072vMts70PwPG1s3IL2HuJ0JsNs5K7GF41MP2rA3EJMHEBMP2PANEJMrIBAPArwr0HCPsxuBM5q7G1MrOHsrsrqNKJCx0v4vwHuBwx2F6vKDIB4Nwrux6DGrML4741C5KLI5KB0N4JAvGtqHw1sFO3q5IvGHqHGFIFMHGFIN2vwBE3s14vq3KN4D6FKHCH2FCvwFsH6PEDCLut2PEBMBIFw3CJE1w503MLq32JCFAtu7Gv0xMDErq727u3MJODG72N25qrE3ODEDG5E1q7sJIJu3MDs32PAF6vs3KJw5u1wFM3u1C5I12NIP4v2t4LwHMLw76ruv2tGP2DuvwJIFA1G5OLw1A3MDALqv4LCH0PwvuLqv412r4PMtuH2DKPGFOPu5630163MBuPuNAHuBw5ML0HEPsvqF2vwJ0xwHAN4JKtEx6D45GN450FCBM7A5w707u3Cx6BCx63srOHsPCr6J6tIxAtIFC1KH45Gr4F07O7ux0NKvI7w32J6PG56vu7q1A3KFqv6xsLCx4rsNsxw10FG1KL0PIFq7qPuNI5q3O1O76BuvED01IPOH2vKF256Bwt47q5GB0NqF4DGJGBAPq5Aru1qDOPurI3wJ4JsL6NKvOtIBsxsFCBuNw10JE1EJAP4JIHAJCFC7EHwv2tuJ6F47CrGH0rEx0vKtMBsxu7sx6rs7wDE1GJs54t4PwJGH2tCrEFCJGJMtAtE30vILOLsLuJArsH2507K7OPqLIvGP656r2PGvwPC5wHqFsLGrM7ErGvuPsr452BuvGJOB6PwFw3qPEPq50J4vw1CH27OvO7uNMB21ANOD6r6BE5GNqLur21MB25sDArO3Ex05A50POF47wt254D2r4rI7uHqP4D2FIBwL6vI1OtANO3OvGHIH6F45sBqH4xO3GrsHuPK7q3q7E7uPCB4HqDExMNEvur2vGHsDCtIt0tuJu54rMv6DEt6vMLKruDAN6PGNuxKtO3GvuD0J4525qPuPAF4BuruPCPOLGFGHuv01M3uPwPqFsvIxE1M3416rqFsvqrMtqHuPM523CN6xqL4H6LwFutAt2F4PAFutItuHwBMxE3OvELMrKPAH4Js3ItA1GvExO3utqFMDCP0LOrI7M7MJEtKF2FE7O7IJ2F21ADODG1urA5w1KPG505q7qBu5GPEt2B0BuN6JItM7OrMxqPwP0r0LCD2tCvK36LO1GDwHO7MxCxM5Ovqv4rE5GJONIv4323qDs7wBwJs32DGNMx2LMtu3K3E5w5qNEv4NMrs3A3wvsP2r2D4r4341KD6LsL47M727sHKN0FG54761MD0vO36twNCtqvKx6JIB4vsrsDKvw1MLsxKLG7svOPGFKrMvGDK5E725MFCHABEJE5I52545MJw16H61qLGHCr2JIr25Ox6NuPKFMDM34tAtqt2tuDM7st4HMPsHsxIrq7GvM36v4tOtI5M5uPCtCxK1IF65ErIH6rAH2Bu1C1M1GFKDuJMPA36v27K5MBONsH4HsJ232HM16DO7IN45412PsFKPGBKJCDqB2LCtwBCPMJwvIrKDEtALw1w367qLqP0HAF0rCLCvG1sF4xMHG74HGPGFw52NwvI1MxAtINGHq1Mv6LA3OD432D4BKBqHI5q5qNAJONGLE1IBMFq3O5wF6H0L2N03OxGxOHCFGtGDGJKN27CtG7sP2B2NKPCBEtCFEDABurM1st6FIB2143uHCFIBuxIJG1E7CPKL6FGvGDwJutGH0tMP0D4PIDwNOxutwx03E7ExG7qFIrIDw7OFw7OLq7qt23GBKFKtO5wJEru505Ov412vIJK3CrqP4r0N67KBG1M7Ex0vOtuHwJ2DMJqt65s36LO5ODsJGL0v6DAB4N2PKtKFw1qtwP41w7uNENMtKFC7wHKJqD6DI3KxM7M3EF4vMtwHq1C10BI1OD0Ds1qNq3stKv6F0Pq70xsF4NCNCB6N4HCH2JwHItIHuPCJMFwJIP2NA7ONGxAPw323CBMxOrAHKNOt2N4tI1M74rO5CvqPELOHGrOr2tELIJqt632BsvI1MHIJ6N0PwBwDqDwPOLI7C7GrqPCtAtIr4JAtqJurKFO1qNEvurqNGFM7KxsL6P0H05OP2PwBO3s3M76N0BGDw3I7ursF2FsLO3utIxKPGJEFO5s5C5O5sFML0vwB0PsvqLAJ2rCJw505Cv4DOrKrG3OHO7uJGNOxutwx01qtuPG52F2xAN6tOBOx4LA3O1wJIJEHIDsLG50HMJAJAJsHE1AJAJA3O1CJEJqLAJqLAJIJAJIJAJC3E3Ev0PEHIDsLG50HMJAJAJsHE1AJAJCJAJAJAJAJAJAJAJAJAJAJAJAJAJAJAJAJAJC1A3O1sJuJAJAJAJAJCJAJCJAJ6NAJAJAJOvMJAJAJAJAJ");
	}
	public DnaStrand[] strands;
	public DnaDomain[] domains;
	public int[] domainComplements;
	
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////

	public void autoWrittenDeSerializeCode(){
		String strands_strTmp = ((StringEntry)readField("strands", new StringEntry(""))).getString();
		if (strands_strTmp.length()>0){
			String[] parts123456 = strands_strTmp.split(",");
			strands = new DnaStrand[parts123456.length];
			for(int qqq = 0; qqq < parts123456.length; qqq++){
				strands[qqq]=new DnaStrand(parts123456[qqq]);
		}}
		String domains_strTmp = ((StringEntry)readField("domains", new StringEntry(""))).getString();
		if (domains_strTmp.length()>0){
			String[] parts123456 = domains_strTmp.split(",");
			domains = new DnaDomain[parts123456.length];
			for(int qqq = 0; qqq < parts123456.length; qqq++){
				domains[qqq]=new DnaDomain(parts123456[qqq]);
		}}
		domainComplements = ((IntArrayEntry)readField("domainComplements", new IntArrayEntry(new int[]{}))).getIntArray();
	}
	public void autoWrittenSerializeCode(){
		writeField("strands", new StringEntry(strands!=null?hashAllToString(strands):""));
		writeField("domains", new StringEntry(domains!=null?hashAllToString(domains):""));
		writeField("domainComplements", new IntArrayEntry(domainComplements));
	}

	
}

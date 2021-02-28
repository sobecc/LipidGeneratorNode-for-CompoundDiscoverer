//-----------------------------------------------------------------------------
// Copyright (c) 2017, Thermo Fisher Scientific
// All rights reserved
//-----------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Diagnostics;
using Thermo.Magellan.BL.Data;
using Thermo.Magellan.BL.Processing;
using Thermo.Magellan.BL.Processing.Interfaces;
using Thermo.Magellan.MassSpec;
using Thermo.Magellan.Utilities;
using Thermo.Metabolism.DataObjects.Constants;

namespace LipidSearch
{
    /// <summary>
    /// The purpose of this node is to generate expected lipids.
    /// </summary>

    #region Node Setup

    [PublisherInformation(Publisher = "Thermo Scientific")]

    [ProcessingNode("9466E17D-5B00-43B1-91CF-93B2C21941CB",
        DisplayName = "Generate Lipids",
        Description = "Generates expected lipids.",
        Category = CDProcessingNodeCategories.UnknownCompounds,
        MainVersion = 1,
        MinorVersion = 1)]

    [ConnectionPoint("Incoming",
        ConnectionDirection = ConnectionDirection.Incoming,
        ConnectionMultiplicity = ConnectionMultiplicity.Single,
        ConnectionMode = ConnectionMode.OnlyImplicitConnectionToAllPossibleParents,
        ConnectionRequirement = ConnectionRequirement.RequiredAtProcessingTime,
        ConnectionDisplayName = ProcessingNodeCategories.DataInput,
        ConnectionDataHandlingType = ConnectionDataHandlingType.InMemory)]
    [ConnectionPointDataContract(
        "Incoming",
        MassSpecDataTypes.SpectrumFiles)]

    [ConnectionPoint("Outgoing",
        ConnectionDirection = ConnectionDirection.Outgoing,
        ConnectionMultiplicity = ConnectionMultiplicity.Multiple,
        ConnectionMode = ConnectionMode.Manual,
        ConnectionRequirement = ConnectionRequirement.Optional,
        ConnectionDataHandlingType = ConnectionDataHandlingType.InMemory)]
    [ConnectionPointDataContract(
        "Outgoing",
        LipidSearchDataTypes.LipidCompounds)]



    [ProcessingNodeConstraints(UsageConstraint = UsageConstraint.Unrestricted)]
    [LicenseFeature(Feature = "CompoundDiscoverer_Base", ShowIfNotAvailable = false)]

    [ProcessingNodeAppearance(ImageLargeSource = "node_icon_scoring_and_annotation_32x32.png")]

    #endregion

    public class LipidGeneratorNode
        : ProcessingNode<IList<LipidCompound>>
        , IResultsSink<IList<SpectrumFile>>
    {
        #region Private Members
        #endregion

        #region Node Parameters
        [StringSelectionParameter(
              Category = "Database Parameters",
              DisplayName = "Class",
              SelectionValues = new string[] { "PC", "PC O", "LPC", "LPC O", "PE", "PE O", "LPE", "LPE O", "PS", "LPS", "PI", "LPI", "PG", "LPG", "PA", "PA O", "LPA", "LPA O", "Cer", "CL", "SM", "TAG", "DAG", "MAG", "HexCer" },
              Description = "Select lipid class. PC phophatidylcholine, PE phosphatidylethanolamine, PS phosphatidylserine, PI phosphatidylethanolamine, PC O- ether phosphatidylcholine, PA phosphatidic acid,Cer ceramide.",
              IsMultiSelect = true,
              DefaultValue = "PC",
              Position = 2)]
        public SimpleSelectionParameter<string> LipidCLass;


        [StringSelectionParameter(
    Category = "Database Parameters",
    DisplayName = "Fatty Acyl",
    SelectionValues = new string[] { "FA 12:0", "FA 12:1", "FA 13:0", "FA 13:1", "FA 14:0", "FA 14:1", "FA 14:2", "FA 15:0", "FA 15:1", "FA 15:2", "FA 16:0", "FA 16:1", "FA 16:2", "FA 16:3", "FA 18:0", "FA 18:1", "FA 18:2", "FA 18:3", "FA 18:4", "FA 19:0", "FA 19:1", "FA 19:2", "FA 19:3", "FA 19:4", "FA 20:0", "FA 20:1", "FA 20:2", "FA 20:3", "FA 20:4", "FA 20:5", "FA 21:0", "FA 21:1", "FA 21:2", "FA 21:3", "FA 21:4", "FA 21:5", "FA 22:0", "FA 22:1", "FA 22:2", "FA 22:3", "FA 22:4", "FA 22:5", "FA 22:6" },
    Description = "Select fatty acids contstituents to search",
    IsMultiSelect = true,
    DefaultValue = "FA 16:0",
    Position = 3)]
        public SimpleSelectionParameter<string> FattyAcyl;
        #endregion

        #region Node Essentials

        /// <summary>
        /// Portion of mass spectra received. (Just to trigger this node to start working.)
        /// </summary>
        public void OnResultsSent(IProcessingNode sender, IList<SpectrumFile> result)
        { }

        /// <summary>
        /// One of the connected nodes sent all of its data.
        /// </summary>
        public override void OnParentNodeFinished(IProcessingNode sender, ResultsArguments eventArgs)
        {
            ArgumentHelper.AssertNotNull(eventArgs, "eventArgs");

            try
            {
                SendAndLogTemporaryMessage("⇒ {0} started...", DisplayName);
                var timer = Stopwatch.StartNew();

                // generate expected lipids
                var lipids = GenerateLipids();

                // send lipids to child nodes
                SendResults(lipids);

                timer.Stop();
                SendAndLogMessage("✓ {0} finished after {1}.", DisplayName, StringHelper.GetDisplayString(timer.Elapsed));
            }
            catch (Exception ex)
            {
                SendAndLogErrorMessage("Error: " + ex.Message);
                throw;
            }

            // inform child nodes that processing has finished
            FireProcessingFinishedEvent(new ResultsArguments(LipidSearchDataTypes.LipidCompounds));
        }

        #endregion

        #region Data Processing

        /// <summary>
        /// Generates expected lipid compounds.
        /// </summary>
        ///
        public IList<LipidCompound> GenerateLipids()
        {
            // init timer
            var timer = Stopwatch.StartNew();

            // init container
            var expectedLipids = new List<LipidCompound>();
            //generates FA building blocks
            #region Make building bloks
            BuildingBlockList NAN = new BuildingBlockList("", "", "0", "", "", "");
            BuildingBlockList FAdb = new BuildingBlockList("FAdb", "FA x:y", "c*x + h*(x-y)*2 + o*2", "12-22", "y=function", "");
            BuildingBlockList SPH = new BuildingBlockList("SPH", "d18:1", "c*18 + h*35 + o*1 + n", "18", "1", "2");

            BuildingBlockList FA = new BuildingBlockList("", "", "", "", "", "");
            FA.RemoveAt(0);

            foreach (var bb in FattyAcyl.Values)
            {
                FA.Add(FAdb.Find(x => x.Name.Contains(bb)));
            }

            var allBB = new Dictionary<string, BuildingBlockList>
            {
                [""] = NAN,
                ["FAdb"] = FAdb,
                ["FA"] = FA,
                ["SPH"] = SPH

            };


            #endregion


            #region Generate List of Candidates
            CompositeCompoundList classPC = new CompositeCompoundList("PC", "PC x:y", "PC r1-r2", "FA ", "", "C8H16O4NP", "FA", "FA", "", "", allBB);
            CompositeCompoundList classPCO = new CompositeCompoundList("PC O", "PC O-x:y", "PC O-r1-r2", "FA ", "", "C8H18O3NP", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPC = new CompositeCompoundList("LPC", "LPC x:y", "LPC r1", "FA ", "", "C8H18O5NP", "FA", "", "", "", allBB);
            CompositeCompoundList classLPCO = new CompositeCompoundList("LPC O", "LPC O-x:y", "LPC O-r1", "FA ", "", "C8H20O4NP", "FA", "", "", "", allBB);

            CompositeCompoundList classPE = new CompositeCompoundList("PE", "PE x:y", "PE r1-r2", "FA ", "", "C5H10O4NP", "FA", "FA", "", "", allBB);
            CompositeCompoundList classPEO = new CompositeCompoundList("PE O", "PE O-x:y", "PE O-r1-r2", "FA ", "", "C5H12O3NP", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPE = new CompositeCompoundList("LPE", "LPE x:y", "LPE r1", "FA ", "", "C5H12O5NP", "FA", "", "", "", allBB);
            CompositeCompoundList classLPEO = new CompositeCompoundList("LPE O", "LPE O-x:y", "LPE O-r1", "FA ", "", "C5H14NO4P", "FA", "", "", "", allBB);

            CompositeCompoundList classPS = new CompositeCompoundList("PS", "PS x:y", "PS r1-r2", "FA ", "", "C6H10NO6P", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPS = new CompositeCompoundList("LPS", "LPS x:y", "LPS r1", "FA ", "", "C6H12NO7P", "FA", "", "", "", allBB);

            CompositeCompoundList classPI = new CompositeCompoundList("PI", "PI x:y", "PI r1-r2", "FA ", "", "C9H15O9P", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPI = new CompositeCompoundList("LPI", "LPI x:y", "LPI r1", "FA ", "", "C9H17O10P", "FA", "", "", "", allBB);

            CompositeCompoundList classPA = new CompositeCompoundList("PA", "PA x:y", "PA r1-r2", "FA ", "", "C3H5O4P", "FA", "FA", "", "", allBB);
            CompositeCompoundList classPAO = new CompositeCompoundList("PA O", "PA O-x:y", "PA O-r1-r2", "FA ", "", "C3H7O3P", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPA = new CompositeCompoundList("LPA", "LPA x:y", "LPA r1", "FA ", "", "C3H7O5P", "FA", "", "", "", allBB);
            CompositeCompoundList classLPAO = new CompositeCompoundList("LPA O", "LPA O-x:y", "LPAO r1", "FA ", "", "C3H9O4P", "FA", "", "", "", allBB);

            CompositeCompoundList classPG = new CompositeCompoundList("PG", "PG x:y", "PG r1-r2", "FA ", "", "C6H11O6", "FA", "FA", "", "", allBB);
            CompositeCompoundList classLPG = new CompositeCompoundList("LPG", "LPG x:y", "LPG r1", "FA ", "", "C6H13O7P", "FA", "", "", "", allBB);

            CompositeCompoundList classCL = new CompositeCompoundList("CL", "CL x:y", "CL r1-r2-r3-r4", "FA ", "", "C9H14O9P2", "FA", "FA", "FA", "FA", allBB);

            CompositeCompoundList classCer = new CompositeCompoundList("Cer", "Cer dx:y", "Cer r1-r2", "SPH", "", "", "SPH", "FA", "", "", allBB);
            CompositeCompoundList classSM = new CompositeCompoundList("SM", "SM dx:y", "SM r1-r2", "SPH", "", "C5H12N2O3P", "SPH", "FA", "", "", allBB);
            CompositeCompoundList classHexCer = new CompositeCompoundList("HexCer", "HexCer dx:y", "HexCer r1-r2", "SPH", "", "C6H10O5", "SPH", "FA", "", "", allBB);

            CompositeCompoundList classTAG = new CompositeCompoundList("TAG", "TAG x:y", "TAG r1-r2-r3", "FA ", "", "C3H2", "FA", "FA", "FA", "", allBB);
            CompositeCompoundList classDAG = new CompositeCompoundList("DAG", "DAG x:y", "DAG r1-r2", "FA ", "", "C3H4O", "FA", "FA", "", "", allBB);
            CompositeCompoundList classMAG = new CompositeCompoundList("MAG", "MAG x:y", "MAG r1", "FA ", "", "C3H4O", "FA", "", "", "", allBB);

            var allClasses = new Dictionary<string, CompositeCompoundList>
            {
                ["PC"] = classPC,
                ["PC O"] = classPCO,
                ["LPC"] = classLPC,
                ["LPC O"] = classLPCO,
                ["PE"] = classPE,
                ["PE O"] = classPEO,
                ["LPE"] = classLPE,
                ["LPE O"] = classLPEO,
                ["PS"] = classPS,
                ["LPS"] = classLPS,
                ["PI"] = classPI,
                ["LPI"] = classLPI,
                ["PA"] = classPA,
                ["PA O"] = classPAO,
                ["LPA"] = classLPA,
                ["LPA O"] = classLPAO,
                ["PG"] = classLPG,
                ["LPG"] = classLPG,
                ["CL"] = classCL,
                ["TAG"] = classTAG,
                ["DAG"] = classDAG,
                ["MAG"] = classMAG,
                ["Cer"] = classCer,
                ["HexCer"] = classHexCer,
                ["SM"] = classSM

            };

            foreach (var key in LipidCLass.Values)
            {
                CompositeCompoundList tempClass = allClasses[key];
                foreach (var lipid in tempClass)
                {
                    if (!expectedLipids.Exists(x => x.Name == lipid.SumComp))
                    {
                        LipidCompound temp = new LipidCompound(lipid.Composition, -1);
                        temp.Name = lipid.SumComp;
                        temp.Class = lipid.Type;
                        expectedLipids.Add(temp);
                    }
                }
            }

            #endregion


            foreach (var species in expectedLipids)
            {
                SendAndLogVerboseMessage(MessageLevel.Info, $"Lipid {species.Name} with {species.ElementalCompositionFormula} on the list of candidates.");
            }
            timer.Stop();
            SendAndLogVerboseMessage(MessageLevel.Info, $"Generating {expectedLipids.Count} expected lipids took {StringHelper.GetDisplayString(timer.Elapsed)}.");

            return expectedLipids;
        }

        #endregion
    }

    #region BuildingBlock Classes
    class BuildingBlock
    {
        public string Type { get; private set; }
        public string Name { get; }
        public string Composition { get; }

        // ???: change x,y,z to an array of int?
        public int X { get; }
        public int Y { get; }
        public int Z { get; }

        public double Mass { get; private set; }

        private BuildingBlock() { }
        public BuildingBlock(string type, string nameTemplate, string compositionTemplate, int x, int y, int z)
        {
            Type = type;
            Name = ExtractInfo.getName(nameTemplate, x, y, z); // make the name (nameTemplate, x, y, z)
            Composition = ExtractInfo.getElemCompFA(compositionTemplate, x, y, z); ; // make the composition (compositionTemplate, x, y, z)
            X = x;
            Y = y;
            Z = z;
            // Mass = ExtractInfo.computeMass(Composition, x, y, z);
        }

        private BuildingBlock(string type, string name, string composition, int x, int y, int z, double mass)
        {
            Type = type;
            Name = name;
            Composition = composition;
            X = x;
            Y = y;
            Z = z;
            Mass = mass;
        }

        override public string ToString() => Name + " with elemental composition " + Composition;

        public static BuildingBlock operator +(BuildingBlock a, BuildingBlock b)
        {
            string type = "type";
            string name = "name";
            string composition = "composition";
            return new BuildingBlock(type, name, composition, a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.Mass + b.Mass);
        }
    }


    class BuildingBlockList : List<BuildingBlock>
    {

        public BuildingBlockList(string tableCompType, string tableCompName, string tableElemComp, string tableX, string tableY, string tableZ)
        {
            List<int> xRanges = ExtractInfo.findIndex(tableX, 0, 0);

            foreach (int x in xRanges)
            {
                List<int> yRanges = ExtractInfo.findIndex(tableY, x, 0);
                foreach (int y in yRanges)
                {
                    List<int> zRanges = ExtractInfo.findIndex(tableZ, x, y);
                    foreach (int z in zRanges)
                    {
                        this.Add(new BuildingBlock(tableCompType, tableCompName, tableElemComp, x, y, z));
                    }

                }
            }
        }

        override public string ToString() => String.Join("\n", values: this);
    }
    #endregion

    #region CompositeCompound Classes
    class CompositeCompound
    {
        public string Type { get; }
        public string SumComp { get; }
        public string MolComp { get; }
        public string Composition { get; }
        public double Mass { get; }

        private CompositeCompound() { }
        public CompositeCompound(string type, string sumCompTemplate, string molCompTemplate, string toDeleteInMolComp, string compositionTemplate, string r1Name, string r2Name, string r3Name, string r4Name, int x, int y, int z)
        {
            Type = type;
            SumComp = ExtractInfo.getName(sumCompTemplate, x, y, z);
            MolComp = getMolComp(molCompTemplate, toDeleteInMolComp, r1Name, r2Name, r3Name, r4Name);
            Composition = ExtractInfo.getElemCompFA(compositionTemplate, x, y, z);

        }

        override public string ToString() => SumComp + " MolComp: " + MolComp + " with elemcomposition " + Composition;

        public string getMolComp(string molCompTemplate, string toDeleteInMolComp, string r1Name, string r2Name, string r3Name, string r4Name)
        {
            string tempMolComp = ExtractInfo.getMolComp(molCompTemplate, r1Name, r2Name, r3Name, r4Name);
            return tempMolComp = tempMolComp.Replace(toDeleteInMolComp, "");
        }

    }
    class CompositeCompoundList : List<CompositeCompound>
    {
        public CompositeCompoundList(string tableCompType, string tableSumCompName, string tableMolCompName, string tableToDeleteInMolComp, string tableConstPart, string tableConstMod, string radical1, string radical2, string radical3, string radical4, Dictionary<string, BuildingBlockList> allBuildingBlocks)
        {

            int limitR1R2; int limitR2R3; int limitR3R4;


            BuildingBlockList R1 = allBuildingBlocks[radical1];
            for (int i = 0; i < R1.Count; i++)
            {

                BuildingBlockList R2 = allBuildingBlocks[radical2];
                if (radical1 == radical2)
                { limitR1R2 = i + 1; }
                else
                {
                    limitR1R2 = R2.Count;
                }
                for (int j = 0; j < limitR1R2; j++)
                {
                    BuildingBlockList R3 = allBuildingBlocks[radical3];
                    if (radical2 == radical3)
                    { limitR2R3 = j + 1; }
                    else
                    {
                        limitR2R3 = R3.Count;
                    }

                    for (int g = 0; g < limitR2R3; g++)
                    {
                        BuildingBlockList R4 = allBuildingBlocks[radical4];
                        if (radical3 == radical4)
                        { limitR3R4 = g + 1; }
                        else
                        {
                            limitR3R4 = R4.Count;
                        }

                        for (int f = 0; f < limitR3R4; f++)
                        {
                            string totalElemComp = R1[i].Composition + "+" + R2[j].Composition + "+" + R3[g].Composition + "+" + R4[f].Composition + "+" + tableConstMod;
                            totalElemComp = totalElemComp.Replace("+0", "");
                            totalElemComp = totalElemComp.Replace("+", "");
                            int totalX = R1[i].X + R2[j].X + R3[g].X + R4[f].X;
                            int totalY = R1[i].Y + R2[j].Y + R3[g].Y + R4[f].Y;
                            int totalZ = R1[i].Z + R2[j].Z + R3[g].Z + R4[f].Z;


                            this.Add(new CompositeCompound(tableCompType, tableSumCompName, tableMolCompName, tableToDeleteInMolComp, totalElemComp, R1[i].Name, R2[j].Name, R3[g].Name, R4[f].Name, totalX, totalY, totalZ));

                        }
                    }
                }
            }
        }

    }
    #endregion

    #region ExtractInfo Helper Class

    static class ExtractInfo
    {
        public static List<int> findIndex(string key, int a, int b)
        {
            List<int> tempLimits = new List<int>();
            string structureInfo = key;


            if (!string.IsNullOrEmpty(structureInfo))
            {

                string[] rangeInfo = structureInfo.Split(new char[] { ',' });

                foreach (string range in rangeInfo)
                {
                    string tmpRange = range;
                    double doubleTempRangeNum;

                    if (double.TryParse(tmpRange, out doubleTempRangeNum))
                    {
                        int intTempRangeNum = (int)doubleTempRangeNum;
                        tempLimits.Add(intTempRangeNum);
                    }
                    else
                    {
                        if (tmpRange.Contains("="))
                        {
                            if (tmpRange.Contains(";"))
                            {
                                int min = (int)getDependendRange(tmpRange.Split(';')[0], a, b);
                                int max = (int)getDependendRange(tmpRange.Split(';')[1], a, b);

                                tempLimits = addIndexToRange(tempLimits, min, max);

                            }
                            else
                            {
                                int min = 0;
                                int max = (int)getDependendRange(tmpRange, a, b);

                                tempLimits = addIndexToRange(tempLimits, min, max);

                            }

                        }
                        else
                        {
                            if (tmpRange.Contains("-"))
                            {
                                int min = int.Parse(tmpRange.Split('-')[0].ToString());
                                int max = int.Parse(tmpRange.Split('-')[1]);

                                tempLimits = addIndexToRange(tempLimits, min, max);
                            }
                        }
                    }
                }
            }
            else

            {
                tempLimits.Add(0);
            }

            return tempLimits;
        }

        public static string getElemComp(string tableElemComp, int indexX, int indexY, int indexZ)
        {
            string ElemCompOut;
            ElemCompOut = tableElemComp.Replace("x", indexX.ToString());
            ElemCompOut = ElemCompOut.Replace("y", indexY.ToString());
            return ElemCompOut;
        }

        public static string getElemCompFA(string tableElemComp, int indexX, int indexY, int indexZ)
        {
            string ElemCompOut;
            int faHydrogens = (indexX - indexY) * 2;
            ElemCompOut = tableElemComp.Replace("(x-y)*2", faHydrogens.ToString());
            ElemCompOut = ElemCompOut.Replace("x", indexX.ToString());
            ElemCompOut = ElemCompOut.Replace("y", indexY.ToString());
            ElemCompOut = ElemCompOut.Replace("z", indexZ.ToString());
            ElemCompOut = ElemCompOut.Replace(" + ", "");
            ElemCompOut = ElemCompOut.Replace("*", "");
            return ElemCompOut.ToUpper();
        }
        public static double getDependendRange(string key, int arg1, int arg2)
        {
            string tmpRange = key;
            double doubleTempValue;

            if (double.TryParse(tmpRange, out doubleTempValue))
            {
                return doubleTempValue;
            }
            else
            {
                Func<int, int, int> func5 = (x, y) => arg1 / 2 - 5;
                return func5.Invoke(arg1, arg2);
            }
        }

        public static List<int> addIndexToRange(List<int> list, int min, int max)
        {
            for (int i = min; i <= max; i++)
            {
                if (!list.Contains(i))
                    list.Add(i);
            }
            return list;
        }

        public static string getName(string tableInput, int x, int y, int z)
        {
            string formattedTableInput = tableInput.Replace("x", "{0}");
            formattedTableInput = formattedTableInput.Replace("y", "{1}");
            formattedTableInput = formattedTableInput.Replace("z", "{2}");
            formattedTableInput = String.Format(formattedTableInput, x, y, z);
            return formattedTableInput;
        }

        public static string getMolComp(string tableInput, string r1, string r2, string r3, string r4)
        {
            string formattedTableInput = tableInput.Replace("r1", "{0}");
            formattedTableInput = formattedTableInput.Replace("r2", "{1}");
            formattedTableInput = formattedTableInput.Replace("r3", "{2}");
            formattedTableInput = formattedTableInput.Replace("r4", "{3}");
            formattedTableInput = String.Format(formattedTableInput, r1, r2, r3, r4);
            return formattedTableInput;
        }
    }
    #endregion
}
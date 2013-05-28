using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using System.Xml;
using System.Xml.Schema;
using System.Xml.Resolvers;
using System.Xml.Serialization;
using System.Xml.XPath;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace OeitTest
{
    /// <summary>
    /// Summary description for SchemaUnitTest
    /// </summary>
    [TestClass]
    public class SchemaUnitTest
    {
        public SchemaUnitTest()
        {
        }

        private TestContext testContextInstance;
        private bool validationErrorsEncountered = false;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        //
        // You can use the following additional attributes as you write your tests:
        //
        // Use ClassInitialize to run code before running the first test in the class
        // [ClassInitialize()]
        // public static void MyClassInitialize(TestContext testContext) { }
        //
        // Use ClassCleanup to run code after all tests in a class have run
        // [ClassCleanup()]
        // public static void MyClassCleanup() { }
        //
        // Use TestInitialize to run code before running each test 
        [TestInitialize]
        public void MyTestInitialize()
        {
            validationErrorsEncountered = false;
        }
        //
        // Use TestCleanup to run code after each test has run
        // [TestCleanup()]
        // public void MyTestCleanup() { }
        //
        #endregion

        [TestMethod]
        public void TestSampleFile1()
        {
            Validate(@"sample1.xml", false);
        }


        private void Validate(string x, bool expected)
        {
            System.IO.TextReader reader = new System.IO.StreamReader(@"oeit.2013.1.xsd");
            XmlSchema schema = XmlSchema.Read(reader, ValidationEventHandler);

            XmlDocument doc = new XmlDocument();
            doc.Schemas.Add(schema);
            doc.Load(x);
            doc.Validate(ValidationEventHandler);

            Assert.AreEqual(expected, validationErrorsEncountered);
        }

        private void ValidationEventHandler(object o, ValidationEventArgs args)
        {
            if (args.Severity == XmlSeverityType.Error || args.Severity == XmlSeverityType.Warning)
            {
                System.Console.Out.WriteLine(args.Severity + ":[" + args.Exception.LineNumber + ":" + args.Exception.LinePosition + "]  " + args.Message);
                validationErrorsEncountered = true;
            }
        }
    }
}

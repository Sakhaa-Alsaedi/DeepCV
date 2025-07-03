#!/usr/bin/env python3
"""
DeepCV PDF Report Generator

This module creates comprehensive PDF reports for DeepCV analysis results,
including all visualizations and summaries from each pipeline step.


"""

import os
from pathlib import Path
import pandas as pd
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
import matplotlib.pyplot as plt
from datetime import datetime

class DeepCVPDFReport:
    """
    Comprehensive PDF report generator for DeepCV analysis
    """
    
    def __init__(self, output_dir="results"):
        """
        Initialize the PDF report generator
        
        Args:
            output_dir (str): Directory containing analysis results
        """
        self.output_dir = Path(output_dir)
        self.viz_dir = self.output_dir / "step_visualizations"
        
        # Setup styles
        self.styles = getSampleStyleSheet()
        self.title_style = ParagraphStyle(
            'CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=24,
            spaceAfter=30,
            alignment=TA_CENTER,
            textColor=colors.darkblue
        )
        self.heading_style = ParagraphStyle(
            'CustomHeading',
            parent=self.styles['Heading2'],
            fontSize=16,
            spaceAfter=12,
            textColor=colors.darkgreen
        )
        self.subheading_style = ParagraphStyle(
            'CustomSubHeading',
            parent=self.styles['Heading3'],
            fontSize=14,
            spaceAfter=8,
            textColor=colors.darkred
        )
    
    def create_comprehensive_report(self, vcf_basename="", step_results=None, plot_paths=None):
        """
        Create comprehensive PDF report with all analysis results
        
        Args:
            vcf_basename (str): Base name for the VCF file
            step_results (dict): Results from each pipeline step
            plot_paths (dict): Paths to generated plots
            
        Returns:
            str: Path to generated PDF report
        """
        print(f"\n📄 Creating Comprehensive PDF Report")
        
        try:
            # Setup PDF document
            report_path = self.output_dir / f"DeepCV_Analysis_Report{vcf_basename}.pdf"
            doc = SimpleDocTemplate(str(report_path), pagesize=A4,
                                  rightMargin=72, leftMargin=72,
                                  topMargin=72, bottomMargin=18)
            
            # Build story (content)
            story = []
            
            # Title page
            story.extend(self._create_title_page(vcf_basename))
            story.append(PageBreak())
            
            # Executive summary
            story.extend(self._create_executive_summary(step_results))
            story.append(PageBreak())
            
            # Step 1: VCF Parsing
            story.extend(self._create_step1_section(plot_paths.get('step1', {}), step_results.get('step1', {})))
            story.append(PageBreak())
            
            # Step 2: SNP Annotation
            story.extend(self._create_step2_section(plot_paths.get('step2', {}), step_results.get('step2', {})))
            story.append(PageBreak())
            
            # Step 3: Risk Assessment
            story.extend(self._create_step3_section(plot_paths.get('step3', {}), step_results.get('step3', {})))
            story.append(PageBreak())
            
            # Step 4: Network Analysis
            story.extend(self._create_step4_section(plot_paths.get('step4', {}), step_results.get('step4', {})))
            story.append(PageBreak())
            
            # Conclusions and recommendations
            story.extend(self._create_conclusions_section(step_results))
            
            # Build PDF
            doc.build(story)
            
            print(f"   PDF report saved: {report_path}")
            return str(report_path)
            
        except Exception as e:
            print(f"   Error creating PDF report: {e}")
            return None
    
    def _create_title_page(self, vcf_basename):
        """Create title page content"""
        story = []
        
        # Main title
        story.append(Paragraph("DeepCV Analysis Report", self.title_style))
        story.append(Spacer(1, 0.5*inch))
        
        # Subtitle
        story.append(Paragraph(f"Comprehensive Genomic Variant Analysis", self.heading_style))
        story.append(Spacer(1, 0.3*inch))
        
        # File information
        if vcf_basename:
            story.append(Paragraph(f"<b>Input File:</b> {vcf_basename}", self.styles['Normal']))
        story.append(Paragraph(f"<b>Analysis Date:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", self.styles['Normal']))
        story.append(Paragraph(f"<b>Generated by:</b> DeepCV Pipeline v1.0", self.styles['Normal']))
        
        story.append(Spacer(1, 1*inch))
        
        # Pipeline overview
        story.append(Paragraph("Pipeline Overview", self.heading_style))
        
        pipeline_text = """
        The DeepCV (Deep Cardiovascular) pipeline is a comprehensive genomic analysis tool 
        that processes VCF files through four main steps:
        
        <b>Step 1:</b> VCF Parsing and Quality Control - Extracts and filters genomic variants
        <b>Step 2:</b> SNP Annotation - Adds functional annotations using multiple databases
        <b>Step 3:</b> Risk Assessment - Calculates disease association scores for genes
        <b>Step 4:</b> Network Analysis - Creates protein-protein interaction networks
        
        This report provides detailed visualizations and analysis results from each step.
        """
        
        story.append(Paragraph(pipeline_text, self.styles['Normal']))
        
        return story
    
    def _create_executive_summary(self, step_results):
        """Create executive summary section"""
        story = []
        
        story.append(Paragraph("Executive Summary", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Key findings
        summary_text = """
        This report presents the results of a comprehensive genomic variant analysis using the DeepCV pipeline. 
        The analysis identified significant genetic variants and their potential associations with various diseases, 
        particularly focusing on cardiovascular and related conditions.
        """
        
        story.append(Paragraph(summary_text, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Key statistics table
        if step_results:
            story.append(Paragraph("Key Statistics", self.heading_style))
            
            # Create summary table
            table_data = [['Analysis Step', 'Key Metrics']]
            
            for step, results in step_results.items():
                if isinstance(results, dict):
                    metrics = ', '.join([f"{k}: {v}" for k, v in list(results.items())[:3]])
                else:
                    metrics = str(results)[:100] + "..." if len(str(results)) > 100 else str(results)
                
                table_data.append([step.replace('_', ' ').title(), metrics])
            
            table = Table(table_data, colWidths=[2*inch, 4*inch])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 12),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            
            story.append(table)
        
        return story
    
    def _create_step1_section(self, plot_paths, step_results):
        """Create Step 1 section content"""
        story = []
        
        story.append(Paragraph("Step 1: VCF Parsing and Quality Control", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Description
        description = """
        The first step of the DeepCV pipeline involves parsing the input VCF file and performing 
        quality control analysis. This includes extracting variant information, calculating quality 
        metrics, and creating visualizations to understand the data distribution.
        """
        
        story.append(Paragraph(description, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Add plots
        for plot_name, plot_path in plot_paths.items():
            if os.path.exists(plot_path):
                story.append(Paragraph(f"{plot_name.replace('_', ' ').title()}", self.subheading_style))
                
                # Add image
                img = Image(plot_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 0.2*inch))
        
        # Key findings
        if step_results:
            story.append(Paragraph("Key Findings", self.heading_style))
            
            findings_text = """
            • Quality score distributions show the overall data quality
            • Manhattan plot reveals variant distribution across chromosomes
            • Variant type analysis identifies the proportion of SNVs, insertions, and deletions
            • Chromosome distribution highlights regions with high variant density
            """
            
            story.append(Paragraph(findings_text, self.styles['Normal']))
        
        return story
    
    def _create_step2_section(self, plot_paths, step_results):
        """Create Step 2 section content"""
        story = []
        
        story.append(Paragraph("Step 2: SNP Annotation", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Description
        description = """
        The second step adds functional annotations to the identified variants using multiple 
        databases and prediction tools. This enriches the variant data with information about 
        gene impacts, consequence types, and functional predictions.
        """
        
        story.append(Paragraph(description, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Add plots
        for plot_name, plot_path in plot_paths.items():
            if os.path.exists(plot_path):
                story.append(Paragraph(f"{plot_name.replace('_', ' ').title()}", self.subheading_style))
                
                # Add image
                img = Image(plot_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 0.2*inch))
        
        # New features added
        features_file = self.viz_dir / "step2_new_features.txt"
        if features_file.exists():
            story.append(Paragraph("New Annotation Features", self.heading_style))
            
            with open(features_file, 'r') as f:
                features_content = f.read()
            
            story.append(Paragraph(features_content.replace('\n', '<br/>'), self.styles['Normal']))
        
        return story
    
    def _create_step3_section(self, plot_paths, step_results):
        """Create Step 3 section content"""
        story = []
        
        story.append(Paragraph("Step 3: Risk Assessment", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Description
        description = """
        The third step performs genetic risk assessment by calculating disease association scores 
        for identified genes. This analysis identifies genes with potential clinical significance 
        and their associations with various diseases and conditions.
        """
        
        story.append(Paragraph(description, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Add plots
        for plot_name, plot_path in plot_paths.items():
            if os.path.exists(plot_path):
                story.append(Paragraph(f"{plot_name.replace('_', ' ').title()}", self.subheading_style))
                
                # Add image
                img = Image(plot_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 0.2*inch))
        
        # Risk assessment summary
        if step_results:
            story.append(Paragraph("Risk Assessment Summary", self.heading_style))
            
            summary_text = """
            • Risk score distributions reveal the range of disease associations
            • Top genes analysis identifies the most clinically relevant genes
            • Disease category analysis shows the types of conditions associated
            • Gene-disease heatmap visualizes the strength of associations
            """
            
            story.append(Paragraph(summary_text, self.styles['Normal']))
        
        return story
    
    def _create_step4_section(self, plot_paths, step_results):
        """Create Step 4 section content"""
        story = []
        
        story.append(Paragraph("Step 4: Network Analysis", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Description
        description = """
        The final step creates protein-protein interaction networks using STRING database 
        to understand the functional relationships between identified genes. This provides 
        insights into biological pathways and potential therapeutic targets.
        """
        
        story.append(Paragraph(description, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Add plots
        for plot_name, plot_path in plot_paths.items():
            if os.path.exists(plot_path):
                story.append(Paragraph(f"{plot_name.replace('_', ' ').title()}", self.subheading_style))
                
                # Add image
                img = Image(plot_path, width=6*inch, height=4*inch)
                story.append(img)
                story.append(Spacer(1, 0.2*inch))
        
        # Network analysis summary
        if step_results:
            story.append(Paragraph("Network Analysis Summary", self.heading_style))
            
            summary_text = """
            • Comprehensive network shows all gene-disease-PPI relationships
            • Disease-gene-PPI network focuses on high-confidence interactions
            • Network statistics provide quantitative measures of connectivity
            • Protein interaction data reveals functional gene modules
            """
            
            story.append(Paragraph(summary_text, self.styles['Normal']))
        
        return story
    
    def _create_conclusions_section(self, step_results):
        """Create conclusions and recommendations section"""
        story = []
        
        story.append(Paragraph("Conclusions and Recommendations", self.title_style))
        story.append(Spacer(1, 0.2*inch))
        
        # Clinical significance
        story.append(Paragraph("Clinical Significance", self.heading_style))
        
        clinical_text = """
        The DeepCV analysis has identified several genes with potential clinical significance 
        based on their disease association scores and protein interaction networks. These 
        findings may be relevant for:
        
        • Personalized medicine approaches
        • Risk stratification for cardiovascular diseases
        • Identification of potential therapeutic targets
        • Understanding of disease mechanisms
        """
        
        story.append(Paragraph(clinical_text, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Recommendations
        story.append(Paragraph("Recommendations", self.heading_style))
        
        recommendations_text = """
        Based on the analysis results, we recommend:
        
        1. <b>Further Validation:</b> Experimental validation of high-scoring gene-disease associations
        2. <b>Clinical Correlation:</b> Correlation with patient phenotypes and clinical outcomes
        3. <b>Pathway Analysis:</b> Detailed pathway enrichment analysis for identified gene sets
        4. <b>Literature Review:</b> Comprehensive literature review for novel associations
        5. <b>Functional Studies:</b> Functional studies for genes with unknown mechanisms
        """
        
        story.append(Paragraph(recommendations_text, self.styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Limitations
        story.append(Paragraph("Limitations", self.heading_style))
        
        limitations_text = """
        This analysis has several limitations that should be considered:
        
        • Results are based on computational predictions and require experimental validation
        • Disease associations are derived from existing databases and may not be complete
        • Protein interaction data may include false positives
        • Clinical significance depends on individual patient context
        """
        
        story.append(Paragraph(limitations_text, self.styles['Normal']))
        
        return story


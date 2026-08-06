import streamlit as st


def show():
    hero_section()
    action_section()
    workflow_section()
    footer_section()


def hero_section():
    st.markdown('''
    <div class = "hero-card">
    <h1>Genova</h1>
    <h2> AI-Powered Genomic Platform </h2>
    <p> 
       Analyze genomic datasets with a simple,
        modern and interactive interface.
    </P>
    </div>   
    ''', unsafe_allow_html=True)


def action_section():
    st.markdown('''
    <div class = "action-container">
    
      <div class = "action-card">
         <h3>Dataset Analyzer</h3>
         <p>parse and validate genomic datasets</p>
      </div>
    
      <div class = "action-card">
          <h3>VCF Analyzer</h3>
          <p>Parse and visualize genomic VCF files.</p>
      </div>
      
      <div class = "action-card">
           <h3>Demo</h3>
           <p>Explore Genova using a sample dataset.</p>
      </div>
    </div>
    ''', unsafe_allow_html=True)
    pass


def workflow_section():
    st.markdown(
        "<h2 class='workflow-title'>How Genova Works</h2>",
        unsafe_allow_html=True
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
          <div class="workflow-card">
              <div class="circle">1</div>
              <h3>Upload Dataset</h3>
              <p>Upload your genomic dataset or a VCF file.</p>
          </div>
          """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
          <div class="workflow-card">
              <div class="circle">2</div>
              <h3>Analyze</h3>
              <p>Genova validates and processes your data.</p>
          </div>
          """, unsafe_allow_html=True)

    with col3:
        st.markdown("""
          <div class="workflow-card">
              <div class="circle">3</div>
              <h3>Visualize Results</h3>
              <p>Explore clear interactive genomic insights.</p>
          </div>
          """, unsafe_allow_html=True)
def footer_section():
    pass

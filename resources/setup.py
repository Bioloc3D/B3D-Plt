from cx_Freeze import setup, Executable
build_exe_options= {
    "include_files": ["logo_B3D.png"],  # Ajoutez vos dépendances
}

setup(
    name="Bioloc3D",
    version="1.0",
    description="Bioloc3D",
    options={"build_exe": build_exe_options},
    executables=[Executable("BD3-Plt.py")],
)
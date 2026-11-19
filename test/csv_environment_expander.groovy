import java.nio.charset.StandardCharsets
import java.nio.file.Files

/**
 * Expands environment variables in text inputs used by NFTest configurations.
 */
class NFTestCsvTools {
    /**
     * Expand $NAME and ${NAME} references using the environment inherited by
     * Nextflow. Undefined variables fail instead of leaving unusable paths.
     */
    String expandEnvironmentVariables(String content) {
        if (content == null) {
            return null
        }

        def variablePattern = ~/\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))/

        return content.replaceAll(variablePattern) { match, bracedName, bareName ->
            def variableName = bracedName ?: bareName
            def variableValue = System.getenv(variableName)

            if (variableValue == null) {
                throw new IllegalArgumentException(
                    "Environment variable ${variableName} referenced by test input is not defined"
                )
            }

            return variableValue
        }
    }

    /**
     * Expand a CSV beneath the Nextflow work directory.
     */
    String prepareInputCsv(Object sourceCsv, String testName, Object launchDir) {
        return prepareInputFile(sourceCsv, testName, launchDir, 'csv')
    }

    /** Expand a project YAML beneath the Nextflow work directory. */
    String prepareInputYaml(Object sourceYaml, String testName, Object launchDir) {
        return prepareInputFile(sourceYaml, testName, launchDir, 'yaml')
    }

    /** Keep expanded inputs in the work tree for inspection after a failed run. */
    String prepareInputFile(Object source, String testName, Object launchDir,
                            String extension) {
        def sourceFile = new File(source.toString())
        if (!sourceFile.isFile()) {
            throw new IllegalArgumentException("Test input file does not exist: ${sourceFile}")
        }

        def workRoot = new File(System.getenv('NXF_WORK') ?: "${launchDir}/work")
        def inputDirectory = new File(workRoot, 'nftest-input')
        Files.createDirectories(inputDirectory.toPath())

        def safeTestName = testName.replaceAll(/[^A-Za-z0-9_.-]/, '-')
        def expandedFile = File.createTempFile("${safeTestName}-", ".${extension}", inputDirectory)
        def sourceContent = Files.readString(sourceFile.toPath(), StandardCharsets.UTF_8)
        def expandedContent = expandEnvironmentVariables(sourceContent)

        Files.writeString(expandedFile.toPath(), expandedContent, StandardCharsets.UTF_8)

        return expandedFile.absolutePath
    }
}

return new NFTestCsvTools()
